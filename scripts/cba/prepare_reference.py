# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Prepare reference network by ensuring all CBA projects are included.

Modify the input network from the SB to get the CBA reference network.

This script adds the missing projects to create
the complete reference network with all projects included.
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
            "prepare_reference",
            planning_horizons="2030",
            run="test-sector-tyndp",
            configfiles=["config/test/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    def normalize_bus_name(bus):
        """Correct instances of UK00/UKNI to GB00/GBNI."""
        return {"UK00": "GB00", "UKNI": "GBNI"}.get(bus, bus)

    def update_or_add_link(n, bus0, bus1, delta, hurdle_costs):
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
                "Updated link %s: p_nom %.2f -> %.2f (delta %.2f)",
                link_name,
                current,
                new_p_nom,
                delta,
            )
            return

        if delta < 0:
            logger.warning(
                "Skipping negative correction %.2f for missing link %s",
                delta,
                link_name,
            )
            return

        if bus0 not in n.buses.index or bus1 not in n.buses.index:
            logger.warning(
                "Skipping link %s: missing buses (%s, %s)",
                link_name,
                bus0,
                bus1,
            )
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
            capital_cost=0.0,
        )
        logger.info("Added missing link %s with p_nom %.2f", link_name, delta)

    # Load simplified network
    n = pypsa.Network(snakemake.input.network)

    # Load project definitions
    transmission_projects = pd.read_csv(snakemake.input.transmission_projects)

    # Get planning horizons from config
    planning_horizons = int(snakemake.wildcards.planning_horizons)

    logger.debug(f"\n{'=' * 80}")
    logger.debug("PREPARING REFERENCE NETWORK")
    logger.debug(f"{'=' * 80}")
    logger.debug(f"Current planning horizon: {planning_horizons}")

    # Hurdle costs: 0.01 €/MWh (p.20, 104 TYNDP 2024 CBA implementation guidelines)
    hurdle_costs = snakemake.params.hurdle_costs

    if planning_horizons in [2035, 2040]:
        corrections = pd.read_csv(snakemake.input.corrections, index_col=0)
        corrections = corrections[
            ["Correction - Summary Direction 1", "Correction - Summary Direction 2"]
        ].copy()
        for border, row in corrections.iterrows():
            if not isinstance(border, str) or "-" not in border:
                continue
            bus0, bus1 = border.split("-", 1)
            update_or_add_link(
                n,
                bus0,
                bus1,
                row["Correction - Summary Direction 1"],
                hurdle_costs,
            )
            update_or_add_link(
                n,
                bus1,
                bus0,
                row["Correction - Summary Direction 2"],
                hurdle_costs,
            )
    else:
        logger.info(
            "Skipping reference corrections for planning horizon %s",
            planning_horizons,
        )

    # Save reference network with all projects
    n.export_to_netcdf(snakemake.output.network)
