# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Prepare reference network by ensuring all CBA projects are included.

Modify the input network from prepare_cba_base to get the CBA reference network.

This script adds the missing projects (from TOOT list) to create
the complete reference network with all projects included.

TODO: Implement project addition logic - currently just passes through the network.
"""

import logging

import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_reference",
            planning_horizons="2030",
            run="test-sector-tyndp",
            configfiles=["config/test/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load base network (already has fixed capacities and hurdle costs)
    n = pypsa.Network(snakemake.input.network)

    # Load project definitions
    transmission_projects = pd.read_csv(snakemake.input.transmission_projects)
    storage_projects = pd.read_csv(snakemake.input.storage_projects)

    # Get planning horizons from config
    planning_horizons = int(snakemake.wildcards.planning_horizons)

    logger.info(f"\n{'=' * 80}")
    logger.info("PREPARING REFERENCE NETWORK")
    logger.info(f"{'=' * 80}")
    logger.info(f"Current planning horizon: {planning_horizons}")

    # TODO: Add missing TOOT projects that are not already in the SB network
    # This ensures that the reference network has all projects included,
    # so that both MSV extraction and TOOT/PINT use the same baseline.
    #
    # For now, pass through the network unchanged (assumes SB already has all projects)

    # Save reference network with all projects
    n.export_to_netcdf(snakemake.output.network)
