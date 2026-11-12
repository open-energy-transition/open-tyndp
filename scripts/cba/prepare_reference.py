# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Prepare the CBA reference network of a planning horizon.

Applies the capacity corrections from the SB to the CBA reference grid
(see `clean_projects.reconcile_reference_grid`) to the simplified SB network.
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
            horizon="2030",
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)

    horizon = snakemake.params.horizon
    corrections = pd.read_csv(snakemake.input.corrections, index_col="link_id").query(
        "planning_horizon == @horizon"
    )["p_nom"]

    missing = corrections.index.difference(n.links.index)
    if not missing.empty:
        logger.warning(
            f"Skipping reference grid corrections for links not in the network: {', '.join(missing)}"
        )
    corrections = corrections.drop(missing)

    logger.info(f"Reference grid corrections for {horizon} (MW):\n{corrections}")
    n.links.loc[corrections.index, "p_nom"] += corrections

    n.export_to_netcdf(snakemake.output.network)
