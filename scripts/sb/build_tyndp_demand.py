# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP Scenario Building demand profiles for Open-TYNDP.

This script processes demand data from TYNDP 2026 (e.g. electricity market,
electricity prosumer, EV charging, hydrogen zones, thermal energy, synthetic
fuels, ...), selecting one of the populated weather scenario columns for demand
profiles.

Weather Scenario Selection
--------------------------

Each planning horizon Excel file contains 30 climate year columns (labeled
``WS001``-``WS030``, ``WS031``-``WS060``, etc.), only 3 of which contain data
for the corresponding planning horizon (the rest are zero-filled placeholders).
Which 3 are populated is dependent on data package and identical across
demand types for a given planning horizon. The availability is recorded in
`AVAILABLE_WEATHER_SCENARIOS`.

`weather_scenarios_tyndp` selects which of them to model, as an ordered
preference list per planning horizon. The first entry that is actually
populated is used; entries that aren't fall back to the first available one
for that horizon (see :py:func:`scripts._helpers.check_weather_scenario`). The
weather scenario is resolved separately for every planning horizon read, since the
``WSxxx`` numbering restarts per horizon.

Data Availability
-----------------

Demand data is only available for 2030, 2035, 2040 and 2050; other planning
horizons are not supported.

Inputs
------

- `data/tyndp/.../2026/Demand`: TYNDP 2026 demand profiles, with one
  subfolder per planning horizon containing one Excel file per demand type
  (e.g. ``ELECTRICITY_MARKET {pyear}.xlsx``, ``Hydrogen_Zone 1_{pyear}.xlsx``),
  each with one sheet per node.

Outputs
-------

- `resources/demand_tyndp_{demand_type}_{planning_horizons}.csv`: Processed
  demand time series for the specified demand type and planning horizon
"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import (
    configure_logging,
    get_snapshots,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

# Weather scenarios (climate year column indices) that contain data in the TYNDP
# 2026 demand files, per planning horizon.
AVAILABLE_WEATHER_SCENARIOS = {
    2030: [3, 21, 29],
    2035: [32, 37, 59],
    2040: [65, 71, 77],
    2050: [91, 92, 106],
}

# Maps the `demand_type` wildcard to the token prefixing the TYNDP 2026 demand
# file names. Also defines which demand types the workflow knows about, and is
# imported by `rules/sb.smk` to constrain the wildcard and collect targets.
DEMAND_TYPE_MAP = {
    "electricity_market": "ELECTRICITY_MARKET",
    "electricity_prosumer": "ELECTRICITY_PROSUMER",
    "electricity_prosumer_btm": "ELECTRICITY_PROSUMER_BEHIND_THE_METER_FIXED_LOAD",
    "ev_market": "EV_FIXED_LOAD_PROFILES_ELECTRICITY_MARKET",
    "ev_prosumer": "EV_FIXED_LOAD_PROFILES_ELECTRICITY_PROSUMER",
    "h2_z1": "Hydrogen_Zone 1",
    "h2_z2": "Hydrogen_Zone 2",
    "synthetic_fuels": "SYNTHETIC_FUELS",
    "thermal_h2": "Thermal_energy_Hydrogen",
    "thermal_ch4": "Thermal_energy_Methane",
}


def get_weather_scenario(weather_scenarios, pyear):
    """
    Select the weather scenario to use for a given planning year.

    Currently selects the first weather year listed for the given planning
    year. This is a placeholder and should be adapted once the full weather
    year implementation is available in SB.

    Parameters
    ----------
    weather_scenarios : dict
        Mapping of planning year to a list of requested weather scenarios,
        e.g. ``{pyear: [weather_scenario, ...]}``.
    pyear : int
        Planning year for which to select the weather scenario.

    Returns
    -------
    int
        The selected weather scenario. Falls back to the first entry in
        ``AVAILABLE_WEATHER_SCENARIOS[pyear]`` if the originally selected
        scenario is not available.

    Notes
    -----
    TODO: currently just selects first weather year, should be adapted when
    we have the full weather year implementation in SB.
    """
    weather_scenario = weather_scenarios[pyear][0]

    if weather_scenario not in AVAILABLE_WEATHER_SCENARIOS[pyear]:
        fallback_scenario = AVAILABLE_WEATHER_SCENARIOS[pyear][0]
        logger.warning(
            f"Weather scenario WS{weather_scenario:03d} not available for "
            f"planning year {pyear}, falling back to WS{fallback_scenario:03d}"
        )
        weather_scenario = fallback_scenario

    return weather_scenario


def check_snapshot_year(year: int, drop_leap_day: bool) -> None:
    """
    Ensure a leap `year` doesn't leave 29 February in `snapshots`.

    TYNDP 2026 demand profiles always span exactly 365 days (no 29 February
    row) — demand is built directly against the target snapshot year (rather
    than a placeholder year later remapped by `align_demand_to_snapshots`).
    A leap `year` is fine as long as `drop_leap_day` strips 29 February from
    `snapshots` (the default); it's only a problem if `year` is a leap year
    and `drop_leap_day` is disabled, since demand can never supply data for
    29 February and it would silently reindex to NaN.

    Raises
    ------
    ValueError
        If `year` is a leap year and `drop_leap_day` is False.
    """
    is_leap_year = pd.Timestamp(year=year, month=1, day=1).is_leap_year
    if is_leap_year and not drop_leap_day:
        raise ValueError(
            f"Snapshot year {year} is a leap year but `enable.drop_leap_day` "
            "is disabled. TYNDP 2026 demand data always spans 365 days (no "
            "29 February). Enable `enable.drop_leap_day` or configure a "
            "non-leap `snapshots` year."
        )


def multiindex_to_datetimeindex(df: pd.DataFrame, year: int) -> pd.DataFrame:
    """Convert demand MultiIndex ('Date', 'Hour') to a DatetimeIndex and return a DataFrame."""

    df_reset = df.reset_index()

    df_reset["datetime"] = pd.to_datetime(
        df_reset["Date"].str.strip(".")
        + f".{year} "
        + (df_reset["Hour"] - 1).astype(str)
        + ":00",
        format="%d.%m.%Y %H:%M",
    )

    # Set as index and drop the old columns
    df_new = df_reset.set_index("datetime").drop(columns=["Date", "Hour"])

    return df_new


def get_file_path(fn: str, pyear: int, demand_type: str) -> Path:
    """
    Construct file path to the demand Excel file for a given planning year and demand type.

    Filenames follow the pattern ``{prefix}[ |_]{pyear}.xlsx``, where `prefix`
    is the raw token `DEMAND_TYPE_MAP` associates with `demand_type`. Neither
    the separator before the year nor the casing of the prefix is consistent
    across demand types, so the file is looked up by matching on the prefix
    instead of a fixed template.

    Parameters
    ----------
    fn : str
        Path to the base directory containing per-year demand data
        subdirectories.
    pyear : int
        Planning year for which to locate the demand file.
    demand_type : str
        Key identifying the demand type, must be present in
        `DEMAND_TYPE_MAP`.

    Returns
    -------
    Path
        Path to the matching demand Excel file. If multiple files match,
        the first one (sorted alphabetically) is returned.

    Raises
    ------
    KeyError
        If `demand_type` is not a known demand type.
    FileNotFoundError
        If no file matches the demand type for `pyear`.
    """
    prefix = DEMAND_TYPE_MAP[demand_type]
    pyear_dir = Path(fn, str(pyear))
    matches = sorted(pyear_dir.glob(f"{prefix}[ _]{pyear}.xlsx"))

    if not matches:
        raise FileNotFoundError(
            f"No demand file found for demand type '{demand_type}' in {pyear_dir}"
        )
    if len(matches) > 1:
        logger.warning(
            f"Multiple demand files match demand type '{demand_type}' in {pyear_dir}: "
            f"{[m.name for m in matches]}. Using {matches[0].name}."
        )

    return matches[0]


def read_demand_excel(demand_fn: str, weather_scenario: int, year: int) -> pd.DataFrame:
    """Read and process demand data from Excel file for a specific weather scenario."""
    ws_code = f"WS{weather_scenario:03d}"
    try:
        data = pd.read_excel(
            demand_fn,
            header=10,
            index_col=[0, 1],
            sheet_name=None,
            usecols=lambda name: name == "Date" or name == "Hour" or name == ws_code,
            engine="calamine",
        )

        demand = pd.concat(data, axis=1).droplevel(1, axis=1)
        if demand.empty:
            raise ValueError(f"No data found for weather scenario {ws_code}.")

        # Build the DatetimeIndex directly against the target snapshot year
        demand = multiindex_to_datetimeindex(demand, year=year)
        # Rename UK in GB
        demand.columns = demand.columns.str.replace("UK", "GB")
        demand.columns.name = "Bus"

    except Exception as e:
        logger.warning(
            f"Failed to read demand from {demand_fn}, weather scenario {ws_code}: "
            f"{type(e).__name__}: {e}"
        )
        demand = pd.DataFrame()

    return demand


def load_demand(
    fn: str,
    pyear: int,
    demand_type: str,
    weather_scenario: int,
    year: int,
) -> pd.DataFrame:
    """
    Load demand data for a planning horizon and demand type.

    Parameters
    ----------
    fn : str
        Path to the base directory containing per-year demand data
        subdirectories.
    pyear : int
        Planning horizon for which to load demand data.
    demand_type : str
        Key identifying the demand type, must be present in
        `DEMAND_TYPE_MAP`.
    weather_scenario : int
        Climate year column index (e.g. 3 for ``WS003``) to load demand data
        for.
    year : int
        Target snapshot year to build the demand DatetimeIndex against.

    Returns
    -------
    pd.DataFrame
        Demand data for the given planning year and demand type.
    """
    demand_fn = get_file_path(fn, pyear, demand_type)
    logger.info(
        f"Processing '{demand_type}' demand for planning horizon {pyear}, "
        f"weather scenario WS{weather_scenario:03d}: reading {demand_fn.name}"
    )

    return read_demand_excel(demand_fn, weather_scenario, year)


# %%
if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_demand",
            planning_horizons="2040",
            demand_type="electricity_market",
            clusters="all",
            run="NT",
            configfiles="config/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    pyear = int(snakemake.wildcards.planning_horizons)
    demand_type = snakemake.wildcards.demand_type
    weather_scenarios = snakemake.params.weather_scenarios
    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )
    fn = snakemake.input.demand

    year = snapshots[0].year
    check_snapshot_year(year, snakemake.params.drop_leap_day)

    weather_scenario = get_weather_scenario(weather_scenarios, pyear)
    demand = load_demand(fn, pyear, demand_type, weather_scenario, year)

    # Export to CSV
    demand.to_csv(snakemake.output.demand, index=True)
