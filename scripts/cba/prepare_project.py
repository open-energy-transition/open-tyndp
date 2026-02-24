# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Prepare a single CBA project network based on the assigned method (TOOT/PINT).

TOOT removes the project from the reference network, PINT adds the project.
Handles multi-border projects, creates links when needed, and validates capacity changes.
"""

import logging

import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config
from scripts.cba._helpers import annuity, get_link_attrs

logger = logging.getLogger(__name__)


def load_method(methods: pd.DataFrame, project_id: int, planning_horizon: int) -> str:
    row = methods[
        (methods["project_id"] == project_id)
        & (methods["planning_horizon"] == planning_horizon)
    ]
    if row.empty:
        raise ValueError(
            f"Missing CBA method for project {project_id} and horizon {planning_horizon}"
        )
    return str(row["method"].iloc[0]).strip().upper()


def apply_toot(n: pypsa.Network, transmission_project: pd.DataFrame) -> None:
    for _, project in transmission_project.iterrows():
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
                "Project removal would result in negative capacity on %s or %s.",
                link_id,
                reverse_link_id,
            )
            raise ValueError("Cannot remove more capacity than exists in the network.")

        if result_capacity == 0:
            n.remove("Link", link_id)
            logger.debug("Removed link %s (capacity reached zero)", link_id)
        else:
            n.links.loc[link_id, "p_nom"] = result_capacity

        if result_capacity_reverse == 0:
            n.remove("Link", reverse_link_id)
            logger.debug("Removed link %s (capacity reached zero)", reverse_link_id)
        else:
            n.links.loc[reverse_link_id, "p_nom"] = result_capacity_reverse


def apply_pint(
    n: pypsa.Network,
    transmission_project: pd.DataFrame,
    hurdle_costs: float,
    annuity_factor: float,
) -> None:
    for _, project in transmission_project.iterrows():
        bus0 = project["bus0"]
        bus1 = project["bus1"]
        link_id = f"{bus0}-{bus1}-DC"
        reverse_link_id = f"{bus1}-{bus0}-DC"

        capacity = project["p_nom 0->1"]
        capacity_reverse = project["p_nom 1->0"]

        if link_id in n.links.index and reverse_link_id in n.links.index:
            n.links.loc[link_id, "p_nom"] += capacity
            n.links.loc[reverse_link_id, "p_nom"] += capacity_reverse
            continue

        for lid, b0, b1, cap in [
            (link_id, bus0, bus1, capacity),
            (reverse_link_id, bus1, bus0, capacity_reverse),
        ]:
            attrs = get_link_attrs(project, cap, annuity_factor)
            n.add(
                "Link",
                lid,
                bus0=b0,
                bus1=b1,
                carrier="DC",
                p_nom=cap,
                marginal_cost=hurdle_costs,
                **attrs,
            )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_project",
            cba_project="t1",
            planning_horizons="2030",
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)
    transmission_projects = pd.read_csv(snakemake.input.transmission_projects)
    methods = pd.read_csv(snakemake.input.methods)

    cba_project = snakemake.wildcards.cba_project
    project_id = int(cba_project[1:])
    planning_horizon = int(snakemake.wildcards.planning_horizons)

    method = load_method(methods, project_id, planning_horizon)
    hurdle_costs = snakemake.params.hurdle_costs

    costs = pd.read_csv(snakemake.input.costs, index_col=0)
    annuity_factor = annuity(
        costs.at["HVDC inverter pair", "lifetime"],
        costs.at["HVDC inverter pair", "discount rate"],
    )

    transmission_project = transmission_projects[
        transmission_projects["project_id"] == project_id
    ]
    assert not transmission_project.empty, (
        f"Transmission project {project_id} not found."
    )

    if method == "TOOT":
        apply_toot(n, transmission_project)
    elif method == "PINT":
        apply_pint(n, transmission_project, hurdle_costs, annuity_factor)
    else:
        raise ValueError(f"Unknown method {method} for project {project_id}")

    logger.info(
        "Saved %s project network for project %s (%s borders)",
        method,
        project_id,
        len(transmission_project),
    )

    n.export_to_netcdf(snakemake.output.network)
