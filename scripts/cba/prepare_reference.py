# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Prepare unified CBA reference network.

The CBA reference network is built from the SB base grid by adding TOOT projects
that are not already in the reference grid. PINT projects are excluded from the
reference (they will be added one at a time for assessment).

Reference = SB base grid + TOOT projects where in_reference=False
"""

import logging

import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config
from scripts.cba._helpers import get_toot_projects, load_method_assignment

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_cba_reference",
            planning_horizons="2030",
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load simplified network (SB base grid)
    n = pypsa.Network(snakemake.input.network)

    # Load project definitions and method assignments
    transmission_projects = pd.read_csv(snakemake.input.transmission_projects)
    method_assignment = load_method_assignment(snakemake.input.method_assignment)

    # Get wildcards
    planning_horizons = int(snakemake.wildcards.planning_horizons)
    scenario = snakemake.config[
        "tyndp_scenario"
    ]  # NT, DE, GA (set by set_scenario_config)

    # Hurdle costs from config
    hurdle_costs = snakemake.params.hurdle_costs

    logger.debug(f"\n{'=' * 80}")
    logger.debug("PREPARING UNIFIED CBA REFERENCE NETWORK")
    logger.debug(f"{'=' * 80}")
    logger.debug(f"Scenario: {scenario}, Planning horizon: {planning_horizons}")

    # Get TOOT projects for this scenario/horizon
    toot_projects = get_toot_projects(
        transmission_projects, planning_horizons, scenario, method_assignment
    )

    logger.info(
        f"TOOT projects for {scenario}{planning_horizons}: {len(toot_projects)}"
    )

    # Track changes
    projects_added = []
    projects_created = []
    projects_in_base = []

    # Add TOOT projects that are NOT already in the reference grid
    in_ref_col = f"in_reference{planning_horizons}"

    for _, project in toot_projects.iterrows():
        in_reference = project[in_ref_col]

        bus0 = project["bus0"]
        bus1 = project["bus1"]
        link_id = f"{bus0}-{bus1}-DC"
        reverse_link_id = f"{bus1}-{bus0}-DC"

        if not in_reference:
            # TOOT project NOT in base grid - ADD it to reference
            capacity = project["p_nom 0->1"]
            capacity_reverse = project["p_nom 1->0"]

            links_exist = link_id in n.links.index and reverse_link_id in n.links.index

            if links_exist:
                # Add capacity to existing links
                original_capacity = n.links.loc[link_id, "p_nom"]
                original_capacity_reverse = n.links.loc[reverse_link_id, "p_nom"]

                n.links.loc[link_id, "p_nom"] += capacity
                n.links.loc[reverse_link_id, "p_nom"] += capacity_reverse

                projects_added.append(project["project_id"])
                logger.debug(
                    f"Project {project['project_id']} ({project['project_name']}) "
                    f"ADDED to existing links:"
                )
                logger.debug(
                    f"    {link_id}: {original_capacity:.0f} → "
                    f"{n.links.loc[link_id, 'p_nom']:.0f} MW (+{capacity:.0f} MW)"
                )
            else:
                # Create new links
                # Get length from project data (from Trans.Investments sheet)
                length = project.get("length", 0)
                if pd.isna(length) or length == 0:
                    logger.warning(
                        f"Project {project['project_id']} ({project['project_name']}) "
                        f"missing length data, using length=0"
                    )
                    length = 0

                # Set underwater_fraction=0 (conservative assumption)
                # TODO: Could derive from Type of Element in Trans.Investments
                underwater_fraction = 0

                n.add(
                    "Link",
                    link_id,
                    bus0=bus0,
                    bus1=bus1,
                    carrier="DC",
                    p_nom=capacity,
                    marginal_cost=hurdle_costs,
                    length=length,
                    underwater_fraction=underwater_fraction,
                )
                n.add(
                    "Link",
                    reverse_link_id,
                    bus0=bus1,
                    bus1=bus0,
                    carrier="DC",
                    p_nom=capacity_reverse,
                    marginal_cost=hurdle_costs,
                    length=length,
                    underwater_fraction=underwater_fraction,
                )

                projects_created.append(project["project_id"])
                logger.debug(
                    f"Project {project['project_id']} ({project['project_name']}) "
                    f"CREATED new links:"
                )
                logger.debug(f"    {link_id}: {capacity:.0f} MW (new)")
                logger.debug(f"    {reverse_link_id}: {capacity_reverse:.0f} MW (new)")
        else:
            # TOOT project already in base grid
            projects_in_base.append(project["project_id"])
            logger.debug(
                f"Project {project['project_id']} ({project['project_name']}) "
                f"already in base network"
            )

    # Summary
    logger.info(f"\n{'=' * 80}")
    logger.info("CBA REFERENCE PREPARATION SUMMARY")
    logger.info(f"{'=' * 80}")
    logger.info(f" Scenario: {scenario}, Horizon: {planning_horizons}")
    logger.info(f" TOOT projects (total): {len(toot_projects['project_id'].unique())}")
    logger.info(f"   - Already in base grid: {len(set(projects_in_base))}")
    logger.info(f"   - Added to existing links: {len(set(projects_added))}")
    logger.info(f"   - Created new links: {len(set(projects_created))}")
    logger.info(f"{'=' * 80}\n")

    # Save reference network
    n.export_to_netcdf(snakemake.output.network)
    logger.info(f"CBA reference network saved ({scenario}{planning_horizons})")
