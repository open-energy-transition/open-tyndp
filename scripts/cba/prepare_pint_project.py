# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Add a single PINT project to the unified CBA reference network.

Creates project networks for PINT methodology by adding one project at a time
to the unified CBA reference network. Handles multi-border projects, creates
new links when needed, and applies hurdle costs to new DC links.
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
            "prepare_pint_project",
            cba_project="t28",
            planning_horizons="2030",
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load CBA reference network
    n = pypsa.Network(snakemake.input.network)
    transmission_projects = pd.read_csv(snakemake.input.transmission_projects)

    cba_project = snakemake.wildcards.cba_project
    project_id = int(cba_project[1:])

    # Hurdle costs for new DC links
    hurdle_costs = snakemake.params.hurdle_costs

    # Get project data
    transmission_project = transmission_projects[
        transmission_projects["project_id"] == project_id
    ]

    assert not transmission_project.empty, (
        f"Transmission project {project_id} not found."
    )

    logger.debug(
        f"PINT: Adding project {project_id} with {len(transmission_project)} border(s)"
    )

    # Add the project (may have multiple borders)
    for _, project in transmission_project.iterrows():
        bus0 = project["bus0"]
        bus1 = project["bus1"]
        link_id = f"{bus0}-{bus1}-DC"
        reverse_link_id = f"{bus1}-{bus0}-DC"

        capacity = project["p_nom 0->1"]
        capacity_reverse = project["p_nom 1->0"]

        # Check if links already exist
        link_exists = link_id in n.links.index
        reverse_link_exists = reverse_link_id in n.links.index

        if link_exists and reverse_link_exists:
            # Links exist - add capacity to existing links
            original_capacity = n.links.loc[link_id, "p_nom"]
            original_capacity_reverse = n.links.loc[reverse_link_id, "p_nom"]

            n.links.loc[link_id, "p_nom"] += capacity
            n.links.loc[reverse_link_id, "p_nom"] += capacity_reverse

            logger.debug(f"Added capacity to existing links for project {project_id}:")
            logger.debug(
                f"  {link_id}: {original_capacity:.0f} → "
                f"{n.links.loc[link_id, 'p_nom']:.0f} MW (+{capacity:.0f} MW)"
            )
            logger.debug(
                f"  {reverse_link_id}: {original_capacity_reverse:.0f} → "
                f"{n.links.loc[reverse_link_id, 'p_nom']:.0f} MW "
                f"(+{capacity_reverse:.0f} MW)"
            )
        else:
            # Create new links
            length = project.get("length", 0)
            if pd.isna(length):
                length = 0

            n.add(
                "Link",
                link_id,
                bus0=bus0,
                bus1=bus1,
                carrier="DC",
                p_nom=capacity,
                marginal_cost=hurdle_costs,
                length=length,
                underwater_fraction=0,
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
                underwater_fraction=0,
            )

            logger.debug(f"Created new links for project {project_id}:")
            logger.debug(f"  {link_id}: {capacity:.0f} MW (new)")
            logger.debug(f"  {reverse_link_id}: {capacity_reverse:.0f} MW (new)")

    logger.info(
        f"PINT project network saved: added project {project_id} "
        f"({len(transmission_project)} border(s))"
    )

    n.export_to_netcdf(snakemake.output.network)
