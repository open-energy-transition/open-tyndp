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
from scripts.cba._helpers import get_link_attrs

logger = logging.getLogger(__name__)

BUS_NAME_MAP = {"UK00": "GB00", "UKNI": "GBNI"}


def normalize_bus_name(bus):
    """Correct instances of UK00/UKNI to GB00/GBNI."""
    return BUS_NAME_MAP.get(bus, bus)


def update_or_add_link(n, bus0, bus1, delta, hurdle_costs, attrs):
    if pd.isna(delta) or delta == 0:
        return

    bus0 = normalize_bus_name(bus0)
    bus1 = normalize_bus_name(bus1)

    link_name = f"{bus0}-{bus1}-DC"
    if link_name in n.links.index:
        current = n.links.at[link_name, "p_nom"]
        new_p_nom = max(0.0, current + delta)
        n.links.at[link_name, "p_nom"] = new_p_nom
        if "p_nom_max" in n.links.columns:
            n.links.at[link_name, "p_nom_max"] = max(
                n.links.at[link_name, "p_nom_max"], new_p_nom
            )
        if "p_nom_min" in n.links.columns:
            n.links.at[link_name, "p_nom_min"] = min(
                n.links.at[link_name, "p_nom_min"], new_p_nom
            )
        logger.info(
            f"Updated link {link_name}: p_nom {current:.2f} -> {new_p_nom:.2f} (delta {delta:.2f})"
        )
        return

    if delta < 0:
        logger.warning(
            f"Skipping negative correction {delta:.2f} for missing link {link_name}"
        )
        return

    if bus0 not in n.buses.index or bus1 not in n.buses.index:
        logger.warning(f"Skipping link {link_name}: missing buses ({bus0}, {bus1})")
        return

    n.add(
        "Link",
        link_name,
        bus0=bus0,
        bus1=bus1,
        carrier="DC",
        p_nom=delta,
        p_nom_max=delta,
        marginal_cost=hurdle_costs,
        **attrs,
    )
    logger.info(f"Added missing link {link_name} with p_nom {delta:.2f}")


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

    # Load simplified network
    n = pypsa.Network(snakemake.input.network)

    costs = pd.read_csv(snakemake.input.costs, index_col=0)

    # Load project definitions
    transmission_projects = pd.read_csv(snakemake.input.transmission_projects)
    storage_projects = pd.read_csv(snakemake.input.storage_projects)

    # Get planning horizons from config
    planning_horizons = int(snakemake.wildcards.planning_horizons)

    logger.debug(f"\n{'=' * 80}")
    logger.debug("PREPARING REFERENCE NETWORK")
    logger.debug(f"{'=' * 80}")
    logger.debug(f"Current planning horizon: {planning_horizons}")

    # Hurdle costs: 0.01 €/MWh (p.20, 104 TYNDP 2024 CBA implementation guidelines)
    hurdle_costs = snakemake.params.hurdle_costs

    if planning_horizons in [2035, 2040] and snakemake.params.patch_sb_with_annexe:
        logger.info(
            f"Reference corrections were already applied for planning horizon {planning_horizons}, skipping"
        )
    elif planning_horizons in [2035, 2040]:
        logger.info("Patching electrical transmission projects with fixes")
        corrections = pd.read_csv(
            snakemake.input.corrections, quotechar="'", index_col=0
        )
        for border, row in corrections.iterrows():
            if not isinstance(border, str) or "-" not in border:
                continue
            bus0, bus1 = border.split("-", 1)

            match = transmission_projects[
                (
                    (transmission_projects["bus0"] == bus0)
                    & (transmission_projects["bus1"] == bus1)
                )
                | (
                    (transmission_projects["bus0"] == bus1)
                    & (transmission_projects["bus1"] == bus0)
                )
            ]
            project = match.iloc[0] if not match.empty else None

            d1 = row["Correction - Summary Direction 1"]
            d2 = row["Correction - Summary Direction 2"]

            attrs = get_link_attrs(project, costs) if project is not None else {}
            update_or_add_link(n, bus0, bus1, d1, hurdle_costs, attrs)
            update_or_add_link(n, bus1, bus0, d2, hurdle_costs, attrs)
    else:
        logger.info(
            f"Skipping reference corrections for planning horizon {planning_horizons}"
        )

    # Save reference network with all projects
    n.export_to_netcdf(snakemake.output.network)
