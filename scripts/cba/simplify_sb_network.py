# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Simplify scenario building network for CBA analysis.

Creates the base CBA network with:
- Fixed optimal capacities from scenario building
- Hurdle costs on DC links
- Extended primary fuel source capacities
- Merge low and high voltage buses

**Inputs**

- Solved network from scenario building workflow

**Outputs**

- `resources/cba/networks/simple_{planning_horizons}.nc`: Simplified network for CBA
"""

import logging

import pandas as pd
import pypsa
from numpy import inf

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def extend_primary_fuel_sources(
    n: pypsa.Network, tyndp_conventional_carriers: list
) -> None:
    """
    Set infinite capacity for primary fuel source generators.

    Rolling horizon with fixed capacities can artificially limit peak fuel production.
    Primary fuel sources have no capital costs, so unlimited capacity ensures
    sufficient supply without affecting the objective function.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify.
    tyndp_conventional_carriers : list
        List of conventional carrier names from TYNDP data, which may include
        fuel sub-types (e.g., 'oil-light', 'oil-heavy'). These are grouped by
        their base fuel type (e.g., 'oil').

    Returns
    -------
    None
        Network is modified in place.
    """
    # split for 'oil-light', 'oil-heavy', 'oil-shale' -> 'oil'
    primary_fuel_carriers = pd.Series(
        [carrier.split("-")[0] for carrier in tyndp_conventional_carriers]
    ).unique()
    mask = n.generators.carrier.str.contains("|".join(primary_fuel_carriers))
    gen_i = n.generators[mask].index
    n.generators.loc[gen_i, "p_nom"] = inf


def move_bus_carrier_and_cleanup(
    n,
    from_carrier="low voltage",
    busmap=None,
):
    """
    Reassign all components attached to buses of `from_carrier` onto their
    corresponding buses (as given by `busmap`), remove the now-empty
    `from_carrier` buses, and drop any lines/links that end up with
    bus0 == bus1 as a result.

    Also drops redundant load-shedding generators sitting on `from_carrier`
    buses, since after reassignment they would duplicate the load-shedding
    generator already present at the target bus.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify in place.
    from_carrier : str, default "low voltage"
        Bus carrier to remove (e.g. "low voltage").
    busmap : pandas.Series, optional
        Mapping of `from_carrier` bus name -> target bus name onto which its
        components should be reassigned. Defaults to
        `n.buses.loc[n.buses.carrier == from_carrier, "location"]`, i.e.
        each bus is mapped onto the bus matching its `location` field.

    Returns
    -------
    n : pypsa.Network
        The modified network (same object, mutated in place).
    """
    if busmap is None:
        busmap = n.buses.loc[n.buses.carrier == from_carrier, "location"]

    # remove load shedding generators on low voltage buses to avoid having two
    # load shedding generators at market node
    load_shedding_gens = n.generators.index[
        (n.generators.bus.map(n.buses.carrier) == from_carrier)
        & (n.generators.carrier == "load")
    ]
    n.remove("Generator", load_shedding_gens)

    # reassign every component attached to the from_carrier buses
    components = sorted(n.branch_components | n.one_port_components)
    for cname in components:
        c = n.c[cname]
        static = c.static
        if static.empty:
            continue
        for port in c.ports:
            col = f"bus{port}"
            if col not in static:
                continue
            mask = static[col].isin(busmap.index)
            if mask.any():
                static.loc[mask, col] = static.loc[mask, col].map(busmap)

    # remove carrier buses
    n.remove("Bus", busmap.index)

    # drop lines/links where bus0 == bus1 after remapping
    for cname in ["Line", "Link"]:
        static = n.c[cname].static
        if static.empty:
            continue
        to_remove = static.index[static.bus0 == static.bus1]
        if len(to_remove):
            n.remove(cname, to_remove)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "simplify_sb_network",
            planning_horizons="2030",
            run="test-sector-tyndp",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the solved network from scenario building
    n = pypsa.Network(snakemake.input.network)

    # Extend primary fuel sources capacity
    tyndp_conventional_carriers = snakemake.params.tyndp_conventional_carriers
    extend_primary_fuel_sources(n, tyndp_conventional_carriers)

    # TODO: in the case of a perfect foresight network we need to extract a single planning horizon here

    # Fix optimal capacities from scenario building
    n.optimize.fix_optimal_capacities()

    # Add hurdle costs to DC links
    # Hurdle costs: 0.01 €/MWh (p.20, 104 TYNDP 2024 CBA implementation guidelines)
    hurdle_costs = snakemake.params.hurdle_costs
    n.links.loc[n.links.carrier == "DC", "marginal_cost"] = hurdle_costs
    logger.info(f"Applied hurdle costs of {hurdle_costs} EUR/MWh to DC links")
    # TODO: for DE/GA add merging of the two H2 zones
    # TODO: for DE/GA add EV electricity consumption from SB as fixed demand

    # merge low voltage to market node
    move_bus_carrier_and_cleanup(n, from_carrier="low voltage")

    # Save base network
    n.export_to_netcdf(snakemake.output.network)
    logger.info(f"Saved CBA base network to {snakemake.output.network}")
