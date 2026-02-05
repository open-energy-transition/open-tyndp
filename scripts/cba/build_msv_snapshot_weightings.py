# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Generate snapshot weightings for MSV extraction temporal aggregation.

Creates a CSV file with resampled snapshot weightings at the configured
MSV extraction resolution (e.g., "24H", "48H").

**Inputs**

- ``resources/cba/networks/reference_{planning_horizons}.nc``: Reference network

**Outputs**

- ``resources/cba/msv_snapshot_weightings_{planning_horizons}.csv``: Snapshot weightings
"""

import logging

import pypsa

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_msv_snapshot_weightings",
            planning_horizons="2030",
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)
    resolution = snakemake.params.msv_resolution

    if not resolution or "h" not in str(resolution).lower():
        # Output empty CSV for non-hourly resolutions (handled differently)
        logger.info("No hourly resolution specified, creating empty weightings file")
        import pandas as pd
        pd.DataFrame().to_csv(snakemake.output.snapshot_weightings)
    else:
        logger.info(f"Resampling snapshot weightings to {resolution}")

        # Resample snapshot weightings to target resolution
        snapshot_weightings = n.snapshot_weightings.resample(resolution).sum()

        # Drop rows with zero weight (gaps in original snapshots)
        snapshot_weightings = snapshot_weightings[snapshot_weightings["objective"] > 0]

        logger.info(
            f"Generated weightings: {len(n.snapshots)} -> {len(snapshot_weightings)} snapshots"
        )
        snapshot_weightings.to_csv(snakemake.output.snapshot_weightings)
