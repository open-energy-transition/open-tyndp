# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Simplify scenario building network for CBA analysis.

Creates the base CBA network with:
- Fixed optimal capacities from scenario building
- Hurdle costs on DC links
- Extended primary fuel source capacities
- Disabled volume limits

**Inputs**

- Solved network from scenario building workflow

**Outputs**

- ``resources/cba/networks/simple_{planning_horizons}.nc``: Simplified network for CBA
"""

import logging

import pandas as pd
import pypsa
from numpy import inf, isfinite

from scripts._helpers import configure_logging, set_scenario_config
from scripts.cba._helpers import summarize_counts

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


def disable_volume_limits(n: pypsa.Network):
    """
    Disable minimum energy production limits (e_sum_min) for generators and links.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify
    """
    for c in n.components[{"Generator", "Link"}]:
        has_e_sum_min = isfinite(c.static.get("e_sum_min", []))
        if has_e_sum_min.any():
            stats = summarize_counts(c.static.loc[has_e_sum_min, "carrier"])
            logger.info(f"Disabling e_sum_min volume limits of:\n{stats}")
            c.static.loc[has_e_sum_min, "e_sum_min"] = -inf


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
    disable_volume_limits(n)

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

    # Save base network
    n.export_to_netcdf(snakemake.output.network)
    logger.info(f"Saved CBA base network to {snakemake.output.network}")
