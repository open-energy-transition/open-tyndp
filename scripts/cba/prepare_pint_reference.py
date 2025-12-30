# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Prepare PINT reference network for CBA analysis.

For PINT (Put IN one at a Time) methodology, the reference network contains
ONLY projects that are already in the reference grid (in_reference=True).
This is essentially the simplified SB network unchanged, since it already
contains only the reference projects.

Projects with in_reference=False are candidates for PINT assessment and will
be added one at a time in prepare_pint_project.py.
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
            "prepare_pint_reference",
            planning_horizons="2030",
            run="test-sector-tyndp",
            configfiles=["config/test/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load simplified network (already contains only in_reference=True projects)
    n = pypsa.Network(snakemake.input.network)

    # Load project definitions for logging/validation
    transmission_projects = pd.read_csv(snakemake.input.transmission_projects)

    # Get planning horizon from wildcard
    planning_horizons = int(snakemake.wildcards.planning_horizons)

    logger.debug(f"\n{'=' * 80}")
    logger.debug("PREPARING PINT REFERENCE NETWORK")
    logger.debug(f"{'=' * 80}")
    logger.debug(f"Current planning horizon: {planning_horizons}")

    # Count projects by reference status
    in_ref_col = f"in_reference{planning_horizons}"
    projects_in_reference = transmission_projects[
        transmission_projects[in_ref_col] == True  # noqa: E712
    ]["project_id"].nunique()
    projects_candidates = transmission_projects[
        transmission_projects[in_ref_col] == False  # noqa: E712
    ]["project_id"].nunique()

    # Summary
    logger.info(f"\n{'=' * 80}")
    logger.info("PINT REFERENCE PREPARATION SUMMARY")
    logger.info(f"{'=' * 80}")
    logger.info(
        f" Projects in reference grid (in_reference=True): {projects_in_reference}"
    )
    logger.info(f" PINT candidate projects (in_reference=False): {projects_candidates}")
    logger.info(" Reference network unchanged (uses SB output directly)")
    logger.info(f"{'=' * 80}\n")

    # Save reference network (unchanged from simplified SB network)
    n.export_to_netcdf(snakemake.output.network)
    logger.info(
        f"PINT reference network saved (horizon {planning_horizons}, "
        f"{projects_in_reference} projects in reference, "
        f"{projects_candidates} candidates for PINT assessment)"
    )
