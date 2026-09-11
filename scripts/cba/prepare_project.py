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
from scripts.cba._helpers import get_link_attrs, get_storage_attrs

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


def load_method(
    methods_fn: str, project_id: int, project_type: str, planning_horizon: int
) -> str:
    """
    Load the method for a specific project and planning horizon.

    Parameters
    ----------
    methods_fn : str
        Path to the file defining the methods.
    project_id : int
        Project reference ID.
    project_type : str
        Either "transmission" or "storage".
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
        & (methods["project_type"] == project_type)
        & (methods["planning_horizon"] == planning_horizon)
    ]
    if row.empty:
        raise ValueError(
            f"Missing CBA method for {project_type} project {project_id} "
            f"and horizon {planning_horizon}"
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


def apply_toot_transmission(
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


def apply_pint_storage(
    n: pypsa.Network,
    storage_project: pd.Series,
    discount_rate: float,
) -> None:
    """
    Add a new CBA storage project as a Bus/Store/Link triple.

    Creates a dedicated storage bus attached to the project's electricity
    bus, a Store sized in MWh, and two Links (charge and discharge) sized in
    MW, using capacities and costs taken from the storage project row.
    """
    project_id = storage_project["project_id"]
    project_name = storage_project["project_name"]
    carrier = storage_project["carrier"]
    ac_bus = storage_project["bus"]
    storage_bus = f"{ac_bus} cba s{project_id} storage"

    attrs = get_storage_attrs(storage_project, discount_rate)

    if carrier not in n.carriers.index:
        n.add("Carrier", carrier)
    n.add("Bus", storage_bus, location=ac_bus, carrier=carrier)
    n.add(
        "Store",
        storage_bus,
        bus=storage_bus,
        carrier=carrier,
        e_nom=attrs["e_nom"],
        e_cyclic=True,
        capital_cost=attrs["capital_cost_per_mwh"],
    )
    n.add(
        "Link",
        f"{storage_bus} charger",
        bus0=ac_bus,
        bus1=storage_bus,
        carrier=carrier,
        p_nom=attrs["p_nom_charge"],
        efficiency=attrs["efficiency"],
    )
    n.add(
        "Link",
        f"{storage_bus} discharger",
        bus0=storage_bus,
        bus1=ac_bus,
        carrier=carrier,
        p_nom=attrs["p_nom_discharge"],
        efficiency=attrs["efficiency"],
    )
    logger.info(
        "Added storage project %s (%s) at bus %s: %.1f MWh, %.1f/%.1f MW charge/discharge",
        project_id,
        project_name,
        ac_bus,
        attrs["e_nom"],
        attrs["p_nom_charge"],
        attrs["p_nom_discharge"],
    )


def prepare_storage_project(
    n: pypsa.Network, snakemake, project_id: int, method: str
) -> None:
    storage_projects = pd.read_csv(snakemake.input.storage_projects)
    storage_project = storage_projects[storage_projects["project_id"] == project_id]
    assert not storage_project.empty, f"Storage project {project_id} not found."

    if method == "toot":
        raise NotImplementedError(
            f"TOOT method not supported for storage project {project_id}: "
            "no matching reference-grid storage component to remove."
        )
    elif method == "pint":
        apply_pint_storage(
            n,
            storage_project.iloc[0],
            snakemake.params.storage_discount_rate,
        )
    else:
        raise ValueError(f"Unknown method {method} for project {project_id}")

    logger.info("Saved %s project network for storage project %s", method, project_id)


def prepare_transmission_project(
    n: pypsa.Network, snakemake, project_id: int, method: str
) -> None:
    transmission_projects = pd.read_csv(snakemake.input.transmission_projects)
    hurdle_costs = snakemake.params.hurdle_costs
    negative_toot_capacity = snakemake.config["cba"].get(
        "negative_toot_capacity", "zero"
    )
    costs = pd.read_csv(snakemake.input.costs, index_col=0)

    transmission_project = transmission_projects[
        transmission_projects["project_id"] == project_id
    ]
    assert not transmission_project.empty, (
        f"Transmission project {project_id} not found."
    )

    if method == "toot":
        apply_toot_transmission(n, transmission_project, negative_toot_capacity)
    elif method == "pint":
        apply_pint_transmission(n, transmission_project, hurdle_costs, costs)
    else:
        raise ValueError(f"Unknown method {method} for project {project_id}")

    logger.info(
        "Saved %s project network for project %s (%s borders)",
        method,
        project_id,
        len(transmission_project),
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

    cba_project = snakemake.wildcards.cba_project
    is_storage = cba_project.startswith("s")
    project_id = int(cba_project[1:])
    planning_horizon = int(snakemake.wildcards.planning_horizons)
    if planning_horizon not in [2030, 2040]:
        logger.warning(
            "CBA methods are only available for 2030 or 2040. Using 2040 for planning horizon %s.",
            snakemake.wildcards.planning_horizons,
        )
        planning_horizon = 2040

    project_type = "storage" if is_storage else "transmission"
    method = load_method(
        snakemake.input.methods, project_id, project_type, planning_horizon
    )

    if is_storage:
        prepare_storage_project(n, snakemake, project_id, method)
    else:
        prepare_transmission_project(n, snakemake, project_id, method)

    n.export_to_netcdf(snakemake.output.network)
