# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Prepare CBA reference network with all TOOT projects included.

Takes the simplified network and ensures all projects that will be evaluated
with TOOT methodology are present. This creates a common baseline for both
MSV extraction and TOOT/PINT evaluation.

TODO: Implement project addition logic - currently passes through unchanged.
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
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
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

    logger.info(f"Current planning horizon: {planning_horizons}")

    # TODO: Add missing TOOT projects that are not already in the SB network
    # This ensures that the reference network has all TOOT projects included,
    # so that both MSV extraction and TOOT/PINT use the same baseline.
    #
    # For now, pass through the network unchanged (assumes SB already has all TOOT projects)

    # Save reference network with all projects
    n.export_to_netcdf(snakemake.output.network)
