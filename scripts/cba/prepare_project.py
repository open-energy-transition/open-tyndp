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


def apply_toot_transmission(
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


def apply_pint_transmission(
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


def apply_pint_storage(
    n: pypsa.Network,
    storage_project: pd.Series,
    discount_rate: float,
    default_lifetime: float,
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

    attrs = get_storage_attrs(storage_project, discount_rate, default_lifetime)

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

    if method == "TOOT":
        raise NotImplementedError(
            f"TOOT method not supported for storage project {project_id}: "
            "no matching reference-grid storage component to remove."
        )
    elif method == "PINT":
        apply_pint_storage(
            n,
            storage_project.iloc[0],
            snakemake.params.storage_discount_rate,
            snakemake.params.storage_default_lifetime,
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

    if method == "TOOT":
        apply_toot_transmission(n, transmission_project, negative_toot_capacity)
    elif method == "PINT":
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
    methods = pd.read_csv(snakemake.input.methods)

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

    method = load_method(methods, project_id, planning_horizon)

    if is_storage:
        prepare_storage_project(n, snakemake, project_id, method)
    else:
        prepare_transmission_project(n, snakemake, project_id, method)

    n.export_to_netcdf(snakemake.output.network)
