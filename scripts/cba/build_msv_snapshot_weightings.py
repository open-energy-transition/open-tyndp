# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Generate snapshot weightings for MSV extraction temporal aggregation.

Produces a CSV with resampled snapshot weightings at the configured
MSV extraction resolution. Follows the same logic as ``time_aggregation.py``
for the supported resolution formats:

- ``false``: No aggregation, outputs empty CSV
- ``"Nsn"``: Representative snapshots (e.g., "2sn"), outputs empty CSV
  (handled directly by ``set_temporal_aggregation``)
- ``"Nh"``: Hourly resampling (e.g., "24H", "48H"), outputs resampled weightings

**Inputs**

- ``resources/cba/networks/reference_{planning_horizons}.nc``: Reference network

**Outputs**

- ``resources/cba/msv_snapshot_weightings_{planning_horizons}.csv``: Snapshot weightings
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
            "build_msv_snapshot_weightings",
            planning_horizons="2030",
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)
    resolution = snakemake.params.msv_resolution
    drop_leap_day = snakemake.params.drop_leap_day

    # No aggregation or representative snapshots — output empty CSV
    # (set_temporal_aggregation handles "Nsn" directly without a file)
    if not resolution or (isinstance(resolution, str) and "sn" in resolution.lower()):
        logger.info("No hourly resampling needed, creating empty weightings file")
        pd.DataFrame().to_csv(snakemake.output.snapshot_weightings)

    # Hourly resampling (e.g., "24H", "48H")
    elif isinstance(resolution, str) and "h" in resolution.lower():
        offset = resolution.lower()
        logger.info(f"Resampling snapshot weightings every {offset}")

        # Resample years separately to handle non-contiguous years
        years = pd.DatetimeIndex(n.snapshots).year.unique()
        snapshot_weightings = []
        for year in years:
            sws_year = n.snapshot_weightings[n.snapshots.year == year]
            sws_year = sws_year.resample(offset).sum()
            snapshot_weightings.append(sws_year)
        snapshot_weightings = pd.concat(snapshot_weightings)

        # Drop rows with zero weight (gaps in original snapshots)
        zeros_i = snapshot_weightings.query("objective == 0").index
        snapshot_weightings.drop(zeros_i, inplace=True)

        # Handle leap days: redistribute weights to March 1st
        swi = snapshot_weightings.index
        leap_days = swi[(swi.month == 2) & (swi.day == 29)]
        if drop_leap_day and not leap_days.empty:
            for year in leap_days.year.unique():
                year_leap_days = leap_days[leap_days.year == year]
                leap_weights = snapshot_weightings.loc[year_leap_days].sum()
                march_first = pd.Timestamp(year, 3, 1, 0, 0, 0)
                snapshot_weightings.loc[march_first] = leap_weights
            snapshot_weightings = snapshot_weightings.drop(leap_days).sort_index()
            logger.info("Dropped leap day(s), redistributed weights to March 1st")

        snapshot_weightings.to_csv(snakemake.output.snapshot_weightings)

    else:
        raise ValueError(
            f"Unsupported MSV extraction resolution: {resolution!r}. "
            "Use false, 'Nsn' (e.g., '2sn'), or 'Nh' (e.g., '24H', '48H')."
        )
