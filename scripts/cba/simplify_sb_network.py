# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Simplify scenario building network for CBA analysis.

Creates the base CBA network with:
- Fixed optimal capacities from scenario building
- Hurdle costs on DC links
- Extended primary fuel source capacities
- Noisy marginal costs removed for StorageUnit and Store components (optional)

**Inputs**

- Solved network from scenario building workflow

**Outputs**

- ``resources/cba/networks/simple_{planning_horizons}.nc``: Simplified network for CBA
"""

import logging

import pandas as pd
import pypsa
from numpy import inf

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def extend_primary_fuel_sources(n: pypsa.Network, tyndp_conventional_carriers: list):
    """
    Set infinite capacity for primary fuel source generators.

    Rolling horizon with fixed capacities can artificially limit peak fuel production.
    Primary fuel sources have no capital costs, so unlimited capacity ensures
    sufficient supply without affecting the objective function.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify
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
    return n


def remove_noisy_marginal_costs(
    n: pypsa.Network,
    threshold: float = 0.02,
    components: list[str] | None = None,
):
    """
    Remove small noisy marginal costs for selected components (Storage, StorageUnit).

    This removes the random noise injected in solve_network.prepare_network
    (approx 0.01 +/- 0.001) by zeroing marginal_cost values below the threshold.
    """
    if components is None:
        components = ["Store", "StorageUnit"]
    for comp_name in components:
        comp = n.components.get(comp_name)
        if comp is None:
            continue
        if "marginal_cost" not in comp.df:
            continue
        mc = comp.df["marginal_cost"]
        cleaned = mc.where(mc.abs() >= threshold, 0.0)
        changed = (cleaned != mc).sum()
        if changed:
            logger.info(
                "Removed noisy marginal_costs for %s: %d entries set to 0 (threshold=%s)",
                comp_name,
                changed,
                threshold,
            )
        comp.df["marginal_cost"] = cleaned
    return n


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

    # Remove noisy marginal costs for Store and StorageUnit components
    if snakemake.params.get("remove_noisy_costs", False):
        remove_noisy_marginal_costs(
            n,
            threshold=snakemake.params.get("noisy_costs_threshold", 0.02),
            components=["Store", "StorageUnit"],
        )

    # Save simplified network
    n.export_to_netcdf(snakemake.output.network)
    logger.info(f"Saved CBA base network to {snakemake.output.network}")
