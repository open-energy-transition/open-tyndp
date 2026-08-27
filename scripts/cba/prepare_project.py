# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Prepare a single CBA project network based on the assigned method (TOOT/PINT).

TOOT removes the project from the reference network, PINT adds the project.
Handles multi-border projects, creates links when needed, and validates capacity changes.
"""

import logging
import random

import numpy as np
import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config
from scripts.cba._helpers import get_link_attrs

logger = logging.getLogger(__name__)


def check_method(method: str) -> str:
    """
    Normalize and validate a given CBA method name.

    Raises
    ------
    ValueError
        If the normalized value is neither "pint" nor "toot".
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


def get_link_capacity_data(n, project, method="toot"):
    """
    Get link IDs and capacities for a DC link project between bus0 and bus1.

    For the TOOT projects, link IDs are looked up directly in `n.links`. For
    the PINT projects, if no matching link exists in the network yet, a
    placeholder ID is constructed instead (e.g. for links to be created).

    Parameters
    ----------
    n : pypsa.Network
        Network containing the links.
    project : pd.Series
        Project data with fields "bus0", "bus1", "p_nom 0->1", "p_nom 1->0".
    method : {"toot", "pint"}, default "toot"
        Lookup strategy. If "pint", missing links fall back to a
        constructed placeholder ID of the form "{bus0}-{bus1}-DC".

    Returns
    -------
    link_id : str
        Forward link (bus0 -> bus1): index of matching links in `n.links`,
        or a placeholder string if method="pint" and none was found.
    reverse_link_id : str
        Reverse link (bus1 -> bus0), same lookup rules as `link_id`.
    capacity : float
        Forward direction capacity (p_nom 0->1).
    capacity_reverse : float
        Reverse direction capacity (p_nom 1->0).
    """
    bus0 = project["bus0"]
    bus1 = project["bus1"]

    link_id = n.links[
        (n.links.bus0 == bus0) & (n.links.bus1 == bus1) & (n.links.carrier == "DC")
    ].index
    reverse_link_id = n.links[
        (n.links.bus0 == bus1) & (n.links.bus1 == bus0) & (n.links.carrier == "DC")
    ].index

    assert len(link_id) <= 1, (
        f"Expected at most one forward link for {bus0}->{bus1}, found {len(link_id)}."
    )
    assert len(reverse_link_id) <= 1, (
        f"Expected at most one reverse link for {bus1}->{bus0}, found {len(reverse_link_id)}."
    )

    if method.lower() == "pint":
        link_id = link_id[0] if not link_id.empty else f"{bus0}-{bus1}-DC"
        reverse_link_id = (
            reverse_link_id[0] if not reverse_link_id.empty else f"{bus1}-{bus0}-DC"
        )
    else:  # TOOT
        if link_id.empty:
            logger.warning(f"TOOT: no forward link found for {bus0} -> {bus1}.")
            link_id = None
        else:
            link_id = link_id[0]
        if reverse_link_id.empty:
            logger.warning(f"TOOT: no reverse link found for {bus1} -> {bus0}.")
            reverse_link_id = None
        else:
            reverse_link_id = reverse_link_id[0]

    capacity = project["p_nom 0->1"]
    capacity_reverse = project["p_nom 1->0"]

    return link_id, reverse_link_id, capacity, capacity_reverse


def apply_toot(
    n: pypsa.Network,
    transmission_project: pd.DataFrame,
    negative_toot_option: str,
) -> None:

    def _apply_toot_capacity(link_id, capacity, project):
        if link_id is None:
            if capacity != 0:
                logger.warning(
                    "Project %s (border: %s) has TOOT capacity of %.0f MW but no matching "
                    "link was found; capacity change skipped.",
                    project["project_id"],
                    project["border"],
                    capacity,
                )
            return
        result_capacity = n.links.loc[link_id, "p_nom"] - capacity
        if result_capacity < 0:
            logger.warning(
                "Applying TOOT for project %s (%s) would create negative capacity: "
                "%s %.0f -> %.0f MW after removing %.0f MW (policy=%s).",
                project["project_id"],
                project["project_name"],
                link_id,
                n.links.loc[link_id, "p_nom"],
                result_capacity,
                capacity,
                negative_toot_option,
            )
            if negative_toot_option == "break":
                raise ValueError(
                    "Cannot remove more capacity than exists in the network."
                )
            if negative_toot_option == "zero":
                result_capacity = max(result_capacity, 0)
            else:
                raise ValueError(
                    f"Unknown cba.negative_toot_option policy: {negative_toot_option}"
                )
        if result_capacity == 0:
            n.remove("Link", link_id)
            logger.debug("Removed link %s (capacity reached zero)", link_id)
        else:
            n.links.loc[link_id, "p_nom"] = result_capacity

    for _, project in transmission_project.iterrows():
        link_id, reverse_link_id, capacity, capacity_reverse = get_link_capacity_data(
            n, project, method="toot"
        )

        _apply_toot_capacity(link_id, capacity, project)
        _apply_toot_capacity(reverse_link_id, capacity_reverse, project)


def apply_pint_transmission(
    n: pypsa.Network,
    transmission_project: pd.DataFrame,
    hurdle_costs: float,
    costs: pd.DataFrame,
) -> None:
    for _, project in transmission_project.iterrows():
        bus0 = project["bus0"]
        bus1 = project["bus1"]

        link_id, reverse_link_id, capacity, capacity_reverse = get_link_capacity_data(
            n, project, method="pint"
        )

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


def _get_generator_values(
    df_static: pd.Series,
    df_dynamic: pd.DataFrame,
    snapshots: pd.Series,
    pypsa_dynamic_attributes: list,
):
    """
    Returns static / dynamic / null value for generator input attributes that can take either a static or time series value

    Parameters
    ----------
    df_static: pd.Series
        Static components of custom generator projects
    df_dynamic: pd.DataFrame
        Dynamic timeseries components of custom generator projects
    snapshots: pd.Series
        pandas DateTime Index
    pypsa_dynamic_attributes: list
        List of PyPSA attributes that can take a static value or series as inputs
    """

    generator_dict = dict()
    for attribute in pypsa_dynamic_attributes:
        if attribute in df_dynamic.columns:
            generator_dict[attribute] = df_dynamic.loc[snapshots, attribute]
        elif attribute in df_static.index:
            generator_dict[attribute] = df_static[attribute]
        else:
            generator_dict[attribute] = np.nan

    return generator_dict


def apply_pint_generator(
    n: pypsa.Network,
    generator_project_static: pd.Series,
    generator_project_dynamic: pd.DataFrame,
) -> None:
    """
    Apply custom generator projects as PINT

    Parameters
    ----------
    n: pypsa.Network
        Network to modify
    generator_project_static: pd.Series
        Static components of custom generator projects
    generator_project_dynamic: pd.DataFrame
        Dynamic components of custom generator projects

    Returns
    -------
    None
    """

    # Generate random hexcode for assigning color to a new carrier
    # Existing color codes are excluded
    def generate_unique_hex(excluded_colors):
        while True:
            # Generate a random 6-digit hex code
            hex_color = f"#{random.randint(0, 0xFFFFFF):06x}"

            # Check if the code is in the exclusion list
            if hex_color not in excluded_colors:
                return hex_color

    # Add generator project to the network
    for _, project in generator_project_static.iterrows():
        # Add carrier to network if new carrier
        if project.carrier not in n.carriers.index:
            n.add(
                "Carrier",
                project.carrier,
                color=generate_unique_hex(
                    n.carriers.color.tolist()
                ),  # Assign a new color to the newly added carrier
            )

        pypsa_dynamic_attributes = (
            n.components["Generator"]
            .defaults.query(
                "type.str.contains('series') and status.str.contains('Input')"
            )
            .index.tolist()
        )
        generator_dict = _get_generator_values(
            project, generator_project_dynamic, n.snapshots, pypsa_dynamic_attributes
        )

        n.add(
            "Generator",
            f"{project.project_name}_{project.project_id}",
            carrier=project.carrier,
            bus=project.bus,
            p_nom=project.p_nom,
            capital_cost=project.capital_cost,
            **generator_dict,
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
    project_id = int(cba_project[1:])
    planning_horizon = int(snakemake.wildcards.planning_horizons)
    methods_fn = snakemake.input.methods
    hurdle_costs = snakemake.params.hurdle_costs
    negative_toot_capacity = snakemake.config["cba"].get(
        "negative_toot_capacity", "zero"
    )

    transmission_projects = pd.read_csv(snakemake.input.transmission_projects)
    generator_projects_static = pd.read_csv(snakemake.input.generator_projects_static)
    generator_projects_dynamic = pd.read_csv(
        snakemake.input.generator_projects_dynamic, header=[0, 1], index_col=0
    )
    costs = pd.read_csv(snakemake.input.costs, index_col=0)
    n = pypsa.Network(snakemake.input.network)

    transmission_project = transmission_projects[
        transmission_projects["project_id"] == project_id
    ]
    generator_project_static = generator_projects_static[
        generator_projects_static["project_id"] == project_id
    ]

    assert not (transmission_project.empty and generator_project_static.empty), (
        f"Transmission or generator project with {project_id} not found."
    )

    if not generator_projects_dynamic.empty:
        generator_project_dynamic = generator_projects_dynamic[str(project_id)]
        generator_project_dynamic.index = pd.to_datetime(
            generator_project_dynamic.index
        )

    if planning_horizon not in [2030, 2040]:
        logger.warning(
            "CBA methods are only available for 2030 or 2040. Using 2040 for planning horizon %s.",
            snakemake.wildcards.planning_horizons,
        )
        planning_horizon = 2040

    method = load_method(methods_fn, project_id, planning_horizon)

    if method == "toot":
        apply_toot(n, transmission_project, negative_toot_capacity)
    elif method == "pint":
        if not transmission_project.empty:
            apply_pint_transmission(n, transmission_project, hurdle_costs, costs)
        elif not generator_project_static.empty:
            apply_pint_generator(n, generator_project_static, generator_project_dynamic)
    else:
        raise ValueError(f"Unknown method {method} for project {project_id}")

    logger.info(
        "Saved %s project network for project %s (%s borders)",
        method,
        project_id,
        len(transmission_project),
    )

    n.export_to_netcdf(snakemake.output.network)
