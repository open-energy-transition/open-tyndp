# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Prepare base network for CBA analysis.

This script creates the common base network used by both:
1. MSV extraction (perfect foresight with cyclicity enabled)
2. Reference network preparation (rolling horizon with cyclicity disabled for seasonal stores)

Operations performed:
- Fix optimal capacities from scenario building
- Add hurdle costs to DC links

**Inputs**

- Solved network from scenario building workflow

**Outputs**

- ``resources/cba/networks/base_{planning_horizons}.nc``: Base network for CBA
"""

import logging

import pypsa

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_cba_base",
            planning_horizons="2030",
            run="test-sector-tyndp",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the solved network from scenario building
    n = pypsa.Network(snakemake.input.network)

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
