# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP Scenario Building demand profiles for Open-TYNDP.

This script processes demand data from TYNDP 2026 (e.g. electricity market,
electricity prosumer, EV charging, hydrogen zones, thermal energy, synthetic
fuels, ...), for one weather scenario.

Weather Scenario Selection
--------------------------

Each planning horizon Excel file contains 30 climate year columns (labeled
``WS001``-``WS030``, ``WS031``-``WS060``, etc.), only 3 of which contain data
for the corresponding planning horizon (the rest are zero-filled placeholders).
Which 3 are populated is dependent on data package and identical across
demand types for a given planning horizon. The availability is recorded in
`AVAILABLE_WEATHER_SCENARIOS`.

Current implementation selects first weather scenario of a planning horizon,
This needs to be revisit once we have implemented the full weather scenario
implementation for SB.

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
  demand time series for the specified demand type and planning horizon.
  Unit depends on `demand_type` (see `DEMAND_TYPE_UNITS`).
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

# Weather scenarios that contain data in the TYNDP 2026 demand files,
# per planning horizon.
AVAILABLE_WEATHER_SCENARIOS = {
    2030: [3, 21, 29],
    2035: [32, 37, 59],
    2040: [65, 71, 77],
    2050: [91, 92, 106],
}

# Maps the `demand_type` wildcard to the TYNDP 2026 demand
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

# Unit of the raw values as they appear in the TYNDP 2026 Excel files 
# confirmed against the official TYNDP 2026 NT+ TimeSeriesDashboard
DEMAND_TYPE_UNITS = {
    "electricity_market": "MW_e",
    "electricity_prosumer": "MW_e",
    "electricity_prosumer_btm": "MW_e",
    "ev_market": "MW_e",
    "ev_prosumer": "MW_e",
    "h2_z1": "MW_H2",
    "h2_z2": "MW_H2",
    "synthetic_fuels": "MW_H2",
    "thermal_h2": "GJ",
    "thermal_ch4": "GJ",
}


def get_weather_scenario(weather_scenarios, pyear):
    """
    Select the weather scenario to use for a given planning year.

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
        Selected weather scenario. Falls back to the first entry in
        ``AVAILABLE_WEATHER_SCENARIOS[pyear]`` if unavailable.

    Notes
    -----
    Currently always picks the first requested weather scenario; should be
    adapted once the full weather year implementation is available in SB.
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

    TYNDP 2026 demand data always spans 365 days, so demand
    built directly against a leap `year` would be missing that day. Therefore,
    `drop_leap_day` needs to be enabled to strip February 29th from `snapshots`.

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
    """
    Convert a ``(Date, Hour)`` MultiIndex to a DatetimeIndex.

    Parameters
    ----------
    df : pd.DataFrame
        Demand data indexed by ``(Date, Hour)``, with ``Date`` as
        ``"DD.MM."`` strings and ``Hour`` as 1-24.
    year : int
        Year to assign to the resulting DatetimeIndex.

    Returns
    -------
    pd.DataFrame
        `df` reindexed with a DatetimeIndex named ``"datetime"``.
    """

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


def deduplicate_corrected_columns(demand: pd.DataFrame) -> pd.DataFrame:
    """
    Prefer a "<node>_corrected" sheet over its plain "<node>" counterpart.

    Some TYNDP raw Excel files carry a leftover "_corrected"-suffixed sheet
    name for a node (e.g. "CZ00_corrected") which is a data-correction artifact.
    If both a plain and a "_corrected" sheet exist for the same node, the plain
    one is dropped and the corrected data is kept; either way, the column is
    renamed back to the plain node code so it lines up with every other bus
    downstream. A warning is logged whenever this happens, since it means
    the raw data needed this override.

    Parameters
    ----------
    demand : pd.DataFrame
        Demand data with one column per bus/node.

    Returns
    -------
    pd.DataFrame
        `demand` with any "_corrected" columns merged into their plain node
        code.
    """
    corrected_cols = [c for c in demand.columns if c.endswith("_corrected")]
    for col in corrected_cols:
        base = col[: -len("_corrected")]
        if base in demand.columns:
            logger.warning(
                f"Found both '{base}' and '{col}' sheets for the same node; "
                f"using the '{col}' (corrected) data for bus '{base}'."
            )
            demand = demand.drop(columns=[base])
        else:
            logger.warning(
                f"Sheet '{col}' found as a corrected variant of node "
                f"'{base}' with no plain '{base}' sheet present; renaming "
                f"it to '{base}'."
            )
        demand = demand.rename(columns={col: base})

    return demand


def drop_zero_demand_columns(demand: pd.DataFrame) -> pd.DataFrame:
    """
    Drop buses with zero total demand, to shrink the resulting output file.

    Some TYNDP nodes can contain zero demand time series for a given demand type
    (e.g. a country with no H2-based heating). Downstream code must treat a missing
    bus as zero demand rather than requiring every bus to be present, since which
    buses are dropped can differ between planning horizons/weather scenarios for
    the same demand type.

    Parameters
    ----------
    demand : pd.DataFrame
        Demand data with one column per bus/node.

    Returns
    -------
    pd.DataFrame
        `demand` without all-zero-sum columns.
    """
    zero_cols = demand.columns[demand.sum(axis=0) == 0]
    if len(zero_cols):
        logger.info(
            f"Dropping {len(zero_cols)} bus(es) with zero total demand: "
            f"{', '.join(zero_cols)}"
        )
        demand = demand.drop(columns=zero_cols)

    return demand


def read_demand_excel(
    demand_fn: str, weather_scenario: int, year: int, demand_type: str
) -> pd.DataFrame:
    """
    Read demand data for one weather scenario from a TYNDP demand Excel file.

    Parameters
    ----------
    demand_fn : str
        Path to the demand Excel file.
    weather_scenario : int
        Climate year column index to read, e.g. 3 for ``WS003``.
    year : int
        Year to assign to the resulting DatetimeIndex.
    demand_type : str
        Key identifying the demand type, must be present in
        `DEMAND_TYPE_UNITS`. Used only to label the output's unit.

    Returns
    -------
    pd.DataFrame
        Demand indexed by DatetimeIndex, one column per bus, with the
        columns' index named ``"Bus [<unit>]"`` (e.g. ``"Bus [MW_e]"``) so
        the unit survives into the output CSV's corner cell. Empty
        DataFrame if reading or parsing fails.
    """
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

        # Build DatetimeIndex from snapshot year
        demand = multiindex_to_datetimeindex(demand, year=year)
        demand = deduplicate_corrected_columns(demand)
        # Rename UK in GB
        demand.columns = demand.columns.str.replace("UK", "GB")
        demand.columns.name = f"Bus [{DEMAND_TYPE_UNITS[demand_type]}]"
        demand = drop_zero_demand_columns(demand)

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
        Climate year column index to load, e.g. 3 for ``WS003``.
    year : int
        Year to assign to the resulting DatetimeIndex.

    Returns
    -------
    pd.DataFrame
        Demand data for the given planning year and demand type.
    """
    demand_fn = get_file_path(fn, pyear, demand_type)
    logger.info(
        f"Processing '{demand_type}' demand ({DEMAND_TYPE_UNITS[demand_type]}) for "
        f"planning horizon {pyear}, weather scenario WS{weather_scenario:03d}: "
        f"reading {demand_fn.name}"
    )

    return read_demand_excel(demand_fn, weather_scenario, year, demand_type)


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
