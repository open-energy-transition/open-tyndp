# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds nuclear availability profiles from the TYNDP 2026 TimeSeries Dashboards.

The hourly nuclear generation of each market node is normalized by its installed capacity
and resampled weekly, retaining the maximum availability factor per week. This accounts for
maintenance schedules and forced outages, which are not given in the PEMMDB.

The profiles are derived from the weather scenario selected for the planning horizon in
``wscenarios_tyndp``.

Outputs
-------
- `resources/nuclear_p_max_pu_{planning_horizons}.csv`: weekly nuclear availability factor per bus.
"""

import logging
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

from scripts._helpers import (
    configure_logging,
    format_bz_names,
    get_wscenario,
    set_scenario_config,
)
from scripts.sb.clean_tyndp_dashboard_benchmark import (
    CAPACITY_SHEET,
    parse_installed_capacity,
    split_unit_from_category,
)

logger = logging.getLogger(__name__)


def read_nuclear_availability(fn: Path) -> pd.DataFrame | None:
    """
    Read the hourly nuclear availability factors from a TYNDP TimeSeries Dashboard file.

    Parameters
    ----------
    fn : Path
        Path to the TimeSeries Dashboard file of one country.

    Returns
    -------
    pd.DataFrame | None
        Hourly nuclear generation normalized by installed capacity, with market nodes as
        columns. None if the country has no nuclear capacity.
    """
    with pd.ExcelFile(fn, engine="calamine") as file:
        capacity = (
            parse_installed_capacity(
                pd.read_excel(file, sheet_name=CAPACITY_SHEET, header=None)
            )
            .query("carrier == 'Nuclear' and sheet == 'E-Market' and value > 0")
            .set_index("bus")
            .value.mul(1e3)  # GW to MW
        )
        if capacity.empty:
            return None
        df = pd.read_excel(file, sheet_name="E-Market", header=None)

    # Row 5 holds the market zone, row 6 the category, hourly values start in row 10
    carrier, _ = split_unit_from_category(df.iloc[6])
    generation = df.loc[10:, carrier == "Nuclear"].astype(float)
    generation.columns = df.loc[5, generation.columns]
    assert generation.columns.is_unique, f"Several nuclear series per node in {fn}"

    return generation.div(capacity)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_nuclear_profiles",
            planning_horizons="2030",
            configfiles="config/config.tyndp.yaml",
            run="NT",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    planning_horizon = int(snakemake.wildcards.planning_horizons)
    wscenario = get_wscenario(snakemake.params.wscenarios, planning_horizon)
    logger.info(
        f"Building nuclear availability profiles for {planning_horizon}, weather scenario WS{wscenario:03d}"
    )

    dashboard_files = sorted(
        Path(
            snakemake.input.dashboard_dir,
            str(planning_horizon),
            f"Weather scenario {wscenario:03d}",
        ).glob("*_TimeSeriesDashboard_*.xlsx")
    )
    with ProcessPoolExecutor(max_workers=snakemake.threads) as executor:
        profiles = pd.concat(
            [
                df
                for df in executor.map(read_nuclear_availability, dashboard_files)
                if df is not None
            ],
            axis=1,
        )

    # Dashboards cover 8760 hours without a year, hence align them to the snapshot year
    year = pd.Timestamp(snakemake.params.snapshots["start"]).year
    sns = pd.date_range(f"{year}-01-01", f"{year}-12-31 23:00", freq="h")
    sns = sns[~((sns.month == 2) & (sns.day == 29))]

    profiles = (
        profiles.set_axis(sns)
        .rename(columns=format_bz_names)
        .rename_axis(index="snapshot", columns=None)
        .resample("W")
        .transform("max")
        .clip(upper=1.0)
        .sort_index(axis=1)
    )
    profiles.to_csv(snakemake.output[0])
