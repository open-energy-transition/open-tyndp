# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Remove a single transmission project from the TOOT reference network.

Creates project networks for TOOT methodology by removing one project at a time.
Handles multi-border projects, removes links when capacity reaches zero, and
validates against negative capacities.
"""

import logging

import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_toot_project", cba_project="t1")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)
    transmission_projects = pd.read_csv(snakemake.input.transmission_projects)

    cba_project = snakemake.wildcards.cba_project
    project_id = int(cba_project[1:])
    method = pd.read_csv(snakemake.input.methods)

    transmission_project = transmission_projects[
        transmission_projects["project_id"] == project_id
    ]

    assert not transmission_project.empty, (
        f"Transmission project {project_id} not found."
    )

    logger.debug(
        f"Project {project_id} has {len(transmission_project)} transmission projects."
    )

    for method, project in transmission_project.iterrows():
        bus0 = project["bus0"]
        bus1 = project["bus1"]
        link_id = f"{bus0}-{bus1}-DC"
        reverse_link_id = f"{bus1}-{bus0}-DC"

        capacity = project["p_nom 0->1"]
        capacity_reverse = project["p_nom 1->0"]

        result_capacity = n.links.loc[link_id, "p_nom"] - capacity
        result_capacity_reverse = (
            n.links.loc[reverse_link_id, "p_nom"] - capacity_reverse
        )

        if result_capacity < 0 or result_capacity_reverse < 0:
            logger.warning(
                f"Project {project_id} removal would result in negative capacity on link {link_id} or {reverse_link_id}."
            )
            raise ValueError("Cannot remove more capacity than exists in the network.")

        if result_capacity == 0:
            n.remove("Link", link_id)
            logger.debug(f"Removed link {link_id} (capacity reached zero)")
        else:
            n.links.loc[link_id, "p_nom"] = result_capacity

        if result_capacity_reverse == 0:
            n.remove("Link", reverse_link_id)
            logger.debug(f"Removed link {reverse_link_id} (capacity reached zero)")
        else:
            n.links.loc[reverse_link_id, "p_nom"] = result_capacity_reverse

    n.export_to_netcdf(snakemake.output.network)
