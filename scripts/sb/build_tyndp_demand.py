# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP Scenario Building demand profiles for Open-TYNDP.

This script processes demand data from TYNDP 2026 (e.g. electricity market,
electricity prosumer, EV charging, hydrogen zones, thermal energy, synthetic
fuels, ...), using the `weather_year` parameter to select one of the available
climate year columns for demand profiles.

Weather Year Selection
-----------------------

Each planning horizon Excel file contains 30 climate year columns (labelled
``WS001``-``WS030``, ``WS031``-``WS060``, etc.), only 3 of which contain data
for the corresponding planning horizon (the rest are empty/zero-filled
placeholders). Which 3 are populated is the same across demand types for a
given planning horizon, and is configured via `load.weather_year_tyndp`
(mapping planning horizon -> list of valid weather years). `weather_year`
selects which of these to use (e.g. `3` for `WS003`); if it isn't valid for
the requested planning horizon, the script falls back to the first entry of
`weather_years` for that horizon (see
:py:func:`scripts._helpers.check_weather_year`).

Data Availability
-----------------

Demand data is available for 2030, 2035, 2040 and 2050. Missing planning
horizons are linearly interpolated between available data points.

Inputs
------

- `data/tyndp_2026_bundle/Demand`: TYNDP 2026 demand profiles, with one
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
    align_demand_to_snapshots,
    check_weather_year,
    configure_logging,
    get_snapshots,
    interpolate_demand,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

# Arbitrary non-leap placeholder year used to build a DatetimeIndex from the
# demand data's day/hour index. Kept fixed across planning years so that
# interpolation between two planning years aligns on matching dates; the
# actual target year is applied later by `align_demand_to_snapshots`.
REFERENCE_YEAR = 2013


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


def get_available_years(fn: str) -> list[int]:
    """Scan the directory to find which planning years are available."""
    return sorted(
        int(folder.name)
        for folder in Path(fn).iterdir()
        if folder.is_dir() and folder.name.isdigit()
    )


def get_file_path(fn: str, pyear: int, demand_type: str) -> Path:
    """
    Construct file path to the demand Excel file for a given planning year and demand type.

    Filenames follow the pattern ``{demand_type}[ |_]{pyear}.xlsx`` (the
    separator before the year isn't consistent across demand types), so the
    file is looked up by matching on the `demand_type` prefix instead of a
    fixed template.
    """
    pyear_dir = Path(fn, str(pyear))
    matches = sorted(pyear_dir.glob(f"{demand_type}*{pyear}.xlsx"))

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


def read_demand_excel(demand_fn: str, weather_year: int) -> pd.DataFrame:
    """Read and process demand data from Excel file for a specific weather year."""
    try:
        ws_code = f"WS{weather_year:03d}"
        data = pd.read_excel(
            demand_fn,
            header=10,
            index_col=[0, 1],
            sheet_name=None,
            usecols=lambda name: name == "Date" or name == "Hour" or name == ws_code,
            engine="calamine",
        )

        demand = pd.concat(data, axis=1).droplevel(1, axis=1)
        # Reindex to match snapshots
        demand = multiindex_to_datetimeindex(demand, year=REFERENCE_YEAR)
        # Rename UK in GB
        demand.columns = demand.columns.str.replace("UK", "GB")
        demand.columns.name = "Bus"

    except Exception as e:
        logger.warning(
            f"Failed to read demand from {demand_fn}, weather year {ws_code}: "
            f"{type(e).__name__}: {e}"
        )
        demand = pd.DataFrame()

    return demand


def load_single_year(
    fn: str,
    pyear: int,
    weather_year: int,
    demand_type: str,
    weather_years: dict[int, list[int]],
) -> pd.DataFrame:
    """Load demand data for a single planning year."""
    demand_fn = get_file_path(fn, pyear, demand_type)

    valid_weather_years = weather_years.get(pyear)
    if valid_weather_years is None:
        logger.warning(
            f"No valid weather years configured for planning year {pyear}. "
            f"Using requested weather year WS{weather_year:03d} as is."
        )
        valid_weather_years = [weather_year]

    weather_year = check_weather_year(weather_year, valid_weather_years)
    return read_demand_excel(demand_fn, weather_year)


def load_demand(
    fn: str,
    pyear: int,
    weather_year: int,
    demand_type: str,
    weather_years: dict[int, list[int]],
) -> pd.DataFrame:
    """
    Load demand data for a specific planning year, weather year and demand type.

    This function retrieves demand data from a file, either by loading the
    exact year if available or by performing linear interpolation between
    available years. The data is filtered for a specific weather year.

    Parameters
    ----------
    fn : str
        Filepath to the demand data directory.
    pyear : int
        Planning year for which to retrieve demand data.
    weather_year : int
        Weather year used to filter the demand data.
    demand_type : str
        Demand type to load (e.g. "ELECTRICITY_MARKET", "Hydrogen_Zone 1").
    weather_years : dict[int, list[int]]
        Mapping of planning horizon to the weather years known to contain
        data for that horizon (see `load.weather_year_tyndp` in the config).

    Returns
    -------
    pd.DataFrame
        DataFrame containing demand data for the specified planning year,
        weather year and demand type.
    """

    available_years = get_available_years(fn)
    logger.info(f"Available years: {available_years}, Target year: {pyear}")

    # If target year exists in data, load it directly
    if pyear in available_years:
        logger.info(f"Year {pyear} found in available data. Loading directly.")
        return load_single_year(fn, pyear, weather_year, demand_type, weather_years)

    # Target year not available, do linear interpolation
    return interpolate_demand(
        available_years=available_years,
        pyear=pyear,
        load_single_year_func=load_single_year,
        fn=fn,
        weather_year=weather_year,
        demand_type=demand_type,
        weather_years=weather_years,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_demand",
            planning_horizons="2040",
            demand_type="ELECTRICITY_MARKET",
            clusters="all",
            configfiles="config/test/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    pyear = int(snakemake.wildcards.planning_horizons)
    demand_type = snakemake.params.demand_type
    weather_year = int(snakemake.params.weather_year)
    weather_years = {
        int(y): list(ws) for y, ws in snakemake.params.weather_years.items()
    }
    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )
    fn = snakemake.input.demand

    logger.info(
        f"Processing '{demand_type}' demand for target year: {pyear}, weather year: WS{weather_year:03d}"
    )
    demand = load_demand(fn, pyear, weather_year, demand_type, weather_years)

    # Reindex demand to fit to snapshots
    demand = align_demand_to_snapshots(demand, snapshots)

    # Export to CSV
    demand.to_csv(snakemake.output.demand, index=True)
