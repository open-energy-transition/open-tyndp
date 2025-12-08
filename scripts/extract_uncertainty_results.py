# SPDX-FileCopyrightText: Contributors to NGV-IEM project
#
# SPDX-License-Identifier: MIT
"""
Extract results from uncertainty scenarios and consolidates them.

Consolidation yields single values for all interconnections, to be used for restricting flows exogenously.
"""

import logging

import pandas as pd
import pypsa

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "extract_uncertainty_results",
            clusters=all,
            sector_opts="",
            planning_horizons=2030,
            configfiles="config/config.tyndp.yaml",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    dispatches = []
    for n_fp in snakemake.input.networks:
        n = pypsa.Network(n_fp)

        # Extract all relevant links (DC links from or to GB00)
        links_s = n.components.links.static
        relevant_links = links_s.loc[
            (links_s["carrier"] == "DC")
            & ((links_s["bus0"] == "GB00") | (links_s["bus1"] == "GB00"))
        ].index

        # Get dispatches for relevant links
        links_d = n.components.links.dynamic
        links_d = links_d["p0"][relevant_links]

        dispatches.append(links_d)

    # Consolidate the results of the different scenarios by averaging
    # TODO: Need to rethink method of aggregation
    combined = pd.concat(dispatches, axis=1)
    consolidated = combined.T.groupby(level=0).mean().T

    # Calculate line limits on a p.u. basis relative to the capacity of each link
    capacities = n.links.loc[consolidated.columns, "p_nom_opt"]
    consolidated = consolidated.div(capacities, axis=1).abs()

    # Set small values to 0
    consolidated = consolidated.where(lambda x: x > 1e-4, 0)

    # Save to CSV
    consolidated.to_csv(snakemake.output["line_limits"])
