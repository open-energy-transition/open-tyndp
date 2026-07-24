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
from scripts.cba._helpers import get_link_attrs

logger = logging.getLogger(__name__)


def check_method(method: str) -> str:
    """
    Normalize and validate the CBA method name.

    If the method is not recognized as either "pint" or "toot", a ValueError is raised.
    """
    method = method.lower().strip()
    if method not in ["pint", "toot"]:
        raise ValueError(f"Method must be 'pint' or 'toot', got: {method}")
    return method


def load_method(methods_fn: str, project_id: int, planning_horizon: int) -> str:
    """
    Load the method for a specific project and planning horizon.

    Parameters
    ----------
    methods_fn : str
        Path to the file defining the methods.
    project_id : int
        Project reference ID.
    planning_horizon : int
        Planning horizon.

    Returns
    -------
    str
        Method to be used to assess a project at a planning horizon.
    """
    methods = pd.read_csv(methods_fn)
    row = methods[
        (methods["project_id"] == project_id)
        & (methods["planning_horizon"] == planning_horizon)
    ]
    if row.empty:
        raise ValueError(
            f"Missing CBA method for project {project_id} and horizon {planning_horizon}"
        )
    return check_method(row["method"].iloc[0])


def apply_toot(
    n: pypsa.Network,
    transmission_project: pd.DataFrame,
    negative_toot_option: str,
) -> None:
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
                "Applying TOOT for project %s (%s) would create negative capacity: "
                "%s %.0f -> %.0f MW after removing %.0f MW, "
                "%s %.0f -> %.0f MW after removing %.0f MW (policy=%s).",
                project["project_id"],
                project["project_name"],
                link_id,
                n.links.loc[link_id, "p_nom"],
                result_capacity,
                capacity,
                reverse_link_id,
                n.links.loc[reverse_link_id, "p_nom"],
                result_capacity_reverse,
                capacity_reverse,
                negative_toot_option,
            )
            if negative_toot_option == "break":
                raise ValueError(
                    "Cannot remove more capacity than exists in the network."
                )
            if negative_toot_option == "zero":
                result_capacity = max(result_capacity, 0)
                result_capacity_reverse = max(result_capacity_reverse, 0)
            else:
                raise ValueError(
                    f"Unknown cba.negative_toot_option policy: {negative_toot_option}"
                )

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
    costs: pd.DataFrame,
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

        attrs = get_link_attrs(project, costs)
        for lid, b0, b1, cap in [
            (link_id, bus0, bus1, capacity),
            (reverse_link_id, bus1, bus0, capacity_reverse),
        ]:
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

    cba_project = snakemake.wildcards.cba_project
    planning_horizon = int(snakemake.wildcards.planning_horizons)
    methods_fn = snakemake.input.methods
    transmission_projects = pd.read_csv(snakemake.input.transmission_projects)
    n = pypsa.Network(snakemake.input.network)
    hurdle_costs = snakemake.params.hurdle_costs
    costs = pd.read_csv(snakemake.input.costs, index_col=0)

    project_id = int(cba_project[1:])
    if planning_horizon not in [2030, 2040]:
        logger.warning(
            "CBA methods are only available for 2030 or 2040. Using 2040 for planning horizon %s.",
            snakemake.wildcards.planning_horizons,
        )
        planning_horizon = 2040

    method = load_method(methods_fn, project_id, planning_horizon)
    negative_toot_capacity = snakemake.config["cba"].get(
        "negative_toot_capacity", "zero"
    )

    transmission_project = transmission_projects[
        transmission_projects["project_id"] == project_id
    ]
    assert not transmission_project.empty, (
        f"Transmission project {project_id} not found."
    )

    if method == "toot":
        apply_toot(n, transmission_project, negative_toot_capacity)
    elif method == "pint":
        apply_pint(n, transmission_project, hurdle_costs, costs)
    else:
        raise ValueError(f"Unknown method {method} for project {project_id}")

    logger.info(
        "Saved %s project network for project %s (%s borders)",
        method,
        project_id,
        len(transmission_project),
    )

    n.export_to_netcdf(snakemake.output.network)
