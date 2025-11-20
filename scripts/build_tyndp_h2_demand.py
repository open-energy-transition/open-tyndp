# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP Scenario Building hydrogen demand profiles for Open-TYNDP.

This script processes hydrogen demand data from TYNDP 2024, using the
``snapshots`` year as the climatic year (``cyear``) for demand profiles.
The data is filtered and interpolated based on the selected scenario
(Distributed Energy, Global Ambition, or National Trends) and planning horizon.

Climatic Year Selection
-----------------------

The ``snapshots`` year determines the climatic year for demand profiles:

- **DE and GA scenarios**: Must use 1995, 2008, or 2009. If ``snapshots``
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

- ``data/tyndp_2024_bundle/Demand Profiles``: TYNDP 2024 hydrogen demand profiles

Outputs
-------

- ``resources/h2_demand_tyndp_{planning_horizons}.csv``: Processed hydrogen
  demand time series for the specified planning horizon
"""

import logging
from pathlib import Path

import pandas as pd
from _helpers import (
    check_cyear,
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
    demand_fn: str, scenario: str, pyear: int, cyear: int, h2_zone: int
) -> pd.DataFrame:
    """Read and process hydrogen demand data from Excel file for a specific year and h2 zone."""
    try:
        data = pd.read_excel(
            demand_fn,
            header=10,
            index_col=[0, 1],
            sheet_name=None,
            usecols=lambda name: name == "Date" or name == "Hour" or name == int(cyear),
        )

        demand = pd.concat(data, axis=1).droplevel(1, axis=1)
        # Reindex to match snapshots
        demand = multiindex_to_datetimeindex(demand, year=cyear)
        # Rename UK in GB
        demand.columns = demand.columns.str.replace("UK", "GB")
        demand.columns.name = "Bus"

    except Exception as e:
        logger.warning(
            f"Failed to read H2 demand for scenario {scenario}, pyear {pyear}, H2 Zone {h2_zone}: "
            f"{type(e).__name__}: {e}"
        )
        demand = pd.DataFrame()

    return demand


def get_file_path(fn: str, scenario: str, pyear: int, h2_zone: int = None) -> Path:
    """Construct file path for given planning year and zone."""

    if scenario == "NT":
        return Path(
            fn,
            scenario,
            "H2 demand profiles",
            f"H2 {pyear}",
            f"{scenario}_{pyear}.xlsx",
        )
    elif scenario in ["DE", "GA"]:
        return Path(
            fn,
            scenario,
            str(pyear),
            f"H2_ZONE_{h2_zone}.xlsx",
        )


def load_single_year(fn: str, scenario: str, pyear: int, cyear: int) -> pd.DataFrame:
    """Load demand data for a single planning year."""
    if scenario == "NT":
        demand_fn = get_file_path(fn, scenario, pyear)
        demand = read_h2_excel(demand_fn, scenario, pyear, cyear, h2_zone=2)
    elif scenario in ["DE", "GA"]:
        demands = {}
        for h2_zone in [1, 2]:
            demand_fn = get_file_path(fn, scenario, pyear, h2_zone)
            demands[h2_zone] = read_h2_excel(
                demand_fn, scenario, pyear, cyear, h2_zone=h2_zone
            )
            demands[h2_zone].columns = [
                f"{col[:2]} H2 Z{h2_zone}" for col in demands[h2_zone].columns
            ]
        demand = pd.concat(demands, axis=1).droplevel(0, axis=1)

    return demand


def load_h2_demand(fn: str, scenario: str, pyear: int, cyear: int) -> pd.DataFrame:
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
    pyear : int
        Planning year for which to retrieve hydrogen demand data.
    cyear : int
        Climatic year used to filter the demand data.

    Returns
    -------
    pd.DataFrame
        DataFrame containing hydrogen demand data for the specified scenario,
        planning year, and climatic year.
    """

    available_years = get_available_years(fn, scenario)
    logger.info(
        f"Scenario {scenario}: Available years: {available_years}, Target year: {pyear}"
    )

    # If target year exists in data, load it directly
    if pyear in available_years:
        logger.info(f"Year {pyear} found in available data. Loading directly.")
        return load_single_year(fn, scenario, pyear, cyear)

    # Target year not available, do linear interpolation
    return interpolate_demand(
        available_years=available_years,
        pyear=pyear,
        load_single_year_func=load_single_year,
        fn=fn,
        scenario=scenario,
        cyear=cyear,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

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
    pyear = int(snakemake.wildcards.planning_horizons)
    snapshots = get_snapshots(snakemake.params.snapshots)
    cyear = snapshots[0].year
    fn = snakemake.input.h2_demand

    # Check if climatic year is valid for scenario
    cyear = check_cyear(cyear, scenario)

    # Load demand with interpolation
    logger.info(
        f"Processing H2 demand for scenario: {scenario}, "
        f"target year: {pyear}, climate year: {cyear}"
    )
    demand = load_h2_demand(fn, scenario, pyear, cyear)

    # Export to CSV
    demand.to_csv(snakemake.output.h2_demand, index=True)
