# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP Scenario Building hydrogen demand profiles for Open-TYNDP.

This script processes hydrogen demand data from TYNDP 2024, using the
`snapshots` year as the climatic year (`weather_scenario`) for demand profiles.
The data is filtered and interpolated based on the selected scenario
(Distributed Energy, Global Ambition, or National Trends) and planning horizon.

Climatic Year Selection
-----------------------

The `snapshots` year determines the climatic year for demand profiles:

- **DE and GA scenarios**: Must use 1995, 2008, or 2009. If `snapshots`
  is not one of these years, 2009 is used as the default (considered most
  representative).
- **NT scenario**: Must be between 1982 and 2019.

Data Availability by Scenario
------------------------------

The input data has different temporal and spatial coverage depending on scenario:

**Distributed Energy (DE) and Global Ambition (GA)**:
  - Hydrogen zone Z2: Available for 2030, 2040, and 2050
  - Hydrogen zone Z1: Available for 2050 (DE and GA) and 2040 (GA only)

**National Trends (NT)**:
  - Available for 2030 and 2040 only
  - No split into hydrogen zones

Processing
----------

Missing years are linearly interpolated between available data points.

Inputs
------

- `data/tyndp_2024_bundle/Demand Profiles`: TYNDP 2024 hydrogen demand profiles

Outputs
-------

- `resources/h2_demand_tyndp_{planning_horizons}.csv`: Processed hydrogen
  demand time series for the specified planning horizon
"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import (
    align_demand_to_snapshots,
    check_weather_scenarios,
    configure_logging,
    get_snapshots,
    interpolate_demand,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def multiindex_to_datetimeindex(df: pd.DataFrame, year: int) -> pd.DataFrame:
    """Convert hydrogen demand MultiIndex ('Date', 'Hour') to a DatetimeIndex and return a DataFrame."""

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


def get_available_years(fn: str, scenario: str) -> list[int]:
    """Scan the directory to find which planning years are available."""
    available_years = []
    scenario_path = Path(fn) / scenario

    if scenario == "NT":
        # Look for folders like "H2 2030", "H2 2040"
        demand_profiles_path = scenario_path / "H2 demand profiles"
        if demand_profiles_path.exists():
            for folder in demand_profiles_path.iterdir():
                if folder.is_dir() and folder.name.startswith("H2 "):
                    year = int(folder.name.split()[-1])
                    available_years.append(year)

    elif scenario in ["DE", "GA"]:
        # Look for year folders directly under scenario
        for folder in scenario_path.iterdir():
            if folder.is_dir() and folder.name.isdigit():
                available_years.append(int(folder.name))

    return sorted(available_years)


def read_h2_excel(
    demand_fn: str, scenario: str, planning_horizon: int, weather_scenario: int, h2_zone: int
) -> pd.DataFrame:
    """Read and process hydrogen demand data from Excel file for a specific year and h2 zone."""
    try:
        data = pd.read_excel(
            demand_fn,
            header=10,
            index_col=[0, 1],
            sheet_name=None,
            usecols=lambda name: name == "Date" or name == "Hour" or name == int(weather_scenario),
        )

        demand = pd.concat(data, axis=1).droplevel(1, axis=1)
        # Reindex to match snapshots
        demand = multiindex_to_datetimeindex(demand, year=weather_scenario)
        # Rename UK in GB
        demand.columns = demand.columns.str.replace("UK", "GB")
        demand.columns.name = "Bus"

    except Exception as e:
        logger.warning(
            f"Failed to read H2 demand for scenario {scenario}, planning_horizon {planning_horizon}, H2 Zone {h2_zone}: "
            f"{type(e).__name__}: {e}"
        )
        demand = pd.DataFrame()

    return demand


def get_file_path(fn: str, scenario: str, planning_horizon: int, h2_zone: int = None) -> Path:
    """
    Construct file path for given planning year and zone.

    Parameters
    ----------
    fn : str
        Base directory path containing scenario subdirectories.
    scenario : str
        Scenario name.
    planning_horizon : int
        Planning year.
    h2_zone : int, optional
        H2 zone identifier required for "DE" and "GA" scenarios. Default is None.

    Returns
    -------
    Path
        Path to the H2 demand Excel file for the given scenario and year.
    """

    if scenario == "NT":
        return Path(
            fn,
            scenario,
            "H2 demand profiles",
            f"H2 {planning_horizon}",
            f"{scenario}_{planning_horizon}.xlsx",
        )
    elif scenario in ["DE", "GA"]:
        return Path(
            fn,
            scenario,
            str(planning_horizon),
            f"H2_ZONE_{h2_zone}.xlsx",
        )


def load_single_year(fn: str, scenario: str, planning_horizon: int, weather_scenario: int) -> pd.DataFrame:
    """Load demand data for a single planning year."""
    if scenario == "NT":
        demand_fn = get_file_path(fn, scenario, planning_horizon)
        demand = read_h2_excel(demand_fn, scenario, planning_horizon, weather_scenario, h2_zone=2)
        demand.columns = [f"{col[:2]} H2" for col in demand.columns]
    elif scenario in ["DE", "GA"]:
        demands = {}
        for h2_zone in [1, 2]:
            demand_fn = get_file_path(fn, scenario, planning_horizon, h2_zone)
            demands[h2_zone] = read_h2_excel(
                demand_fn, scenario, planning_horizon, weather_scenario, h2_zone=h2_zone
            )
            demands[h2_zone].columns = [
                f"{col[:2]} H2 Z{h2_zone}" for col in demands[h2_zone].columns
            ]
        demand = pd.concat(demands, axis=1).droplevel(0, axis=1)

    return demand


def load_h2_demand(fn: str, scenario: str, planning_horizon: int, weather_scenario: int) -> pd.DataFrame:
    """
    Load hydrogen demand data for a specific scenario, climate year, planning year.

    This function retrieves hydrogen demand data from a file, either by loading
    the exact year if available or by performing linear interpolation between
    available years. The data is filtered for a specific climatic year.

    Parameters
    ----------
    fn : str
        Filepath to the hydrogen demand data file.
    scenario : str
        Name of the scenario to load.
    planning_horizon : int
        Planning year for which to retrieve hydrogen demand data.
    weather_scenario : int
        Climatic year used to filter the demand data.

    Returns
    -------
    pd.DataFrame
        DataFrame containing hydrogen demand data for the specified scenario,
        planning year, and climatic year.
    """

    available_years = get_available_years(fn, scenario)
    logger.info(
        f"Scenario {scenario}: Available years: {available_years}, Target year: {planning_horizon}"
    )

    # If target year exists in data, load it directly
    if planning_horizon in available_years:
        logger.info(f"Year {planning_horizon} found in available data. Loading directly.")
        return load_single_year(fn, scenario, planning_horizon, weather_scenario)

    # Target year not available, do linear interpolation
    return interpolate_demand(
        available_years=available_years,
        planning_horizon=planning_horizon,
        load_single_year_func=load_single_year,
        fn=fn,
        scenario=scenario,
        weather_scenario=weather_scenario,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_h2_demand",
            planning_horizons="2040",
            clusters="all",
            configfiles="config/test/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    scenario = snakemake.params["scenario"]
    planning_horizon = int(snakemake.wildcards.planning_horizons)
    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )
    weather_scenario = snapshots[0].year
    fn = snakemake.input.h2_demand

    if scenario not in ["DE", "GA", "NT"]:
        demand = pd.Series()

    else:
        # Check if climatic year is valid for scenario
        weather_scenario = check_weather_scenarios(weather_scenario, scenario)

        # Load demand with interpolation
        logger.info(
            f"Processing H2 demand for scenario: {scenario}, "
            f"target year: {planning_horizon}, climate year: {weather_scenario}"
        )
        demand = load_h2_demand(fn, scenario, planning_horizon, weather_scenario)

        # Reindex demand to fit to snapshots
        demand = align_demand_to_snapshots(demand, snapshots)

    # Export to CSV
    demand.to_csv(snakemake.output.h2_demand, index=True)
