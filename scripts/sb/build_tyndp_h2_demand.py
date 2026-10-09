# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP Scenario Building hydrogen demand profiles for Open-TYNDP.

This script processes hydrogen demand data from TYNDP 2024, using the
`snapshots` year as the climatic year (`wscenario`) for demand profiles.
The data is filtered and interpolated based on the planning horizon.

Climatic Year Selection
-----------------------

The `snapshots` year determines the climatic year for demand profiles. It must
be 1995, 2008, or 2009. Otherwise, 2009 is used as the default (considered most
representative).

Data Availability
-----------------

The input data covers the National Trends (NT) scenario:
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
    check_wscenario,
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


def get_available_years(fn: str) -> list[int]:
    """Scan the directory to find which planning years are available."""
    available_years = []

    # Look for folders like "H2 2030", "H2 2040"
    demand_profiles_path = Path(fn) / "NT" / "H2 demand profiles"
    if demand_profiles_path.exists():
        for folder in demand_profiles_path.iterdir():
            if folder.is_dir() and folder.name.startswith("H2 "):
                year = int(folder.name.split()[-1])
                available_years.append(year)

    return sorted(available_years)


def read_h2_excel(
    demand_fn: str,
    planning_horizon: int,
    wscenario: int,
) -> pd.DataFrame:
    """Read and process hydrogen demand data from Excel file for a specific year."""
    try:
        data = pd.read_excel(
            demand_fn,
            header=10,
            index_col=[0, 1],
            sheet_name=None,
            usecols=lambda name: (
                name == "Date" or name == "Hour" or name == int(wscenario)
            ),
        )

        demand = pd.concat(data, axis=1).droplevel(1, axis=1)
        # Reindex to match snapshots
        demand = multiindex_to_datetimeindex(demand, year=wscenario)
        # Rename UK in GB
        demand.columns = demand.columns.str.replace("UK", "GB")
        demand.columns.name = "Bus"

    except Exception as e:
        logger.warning(
            f"Failed to read H2 demand for planning_horizon {planning_horizon}: "
            f"{type(e).__name__}: {e}"
        )
        demand = pd.DataFrame()

    return demand


def load_single_year(fn: str, planning_horizon: int, wscenario: int) -> pd.DataFrame:
    """Load demand data for a single planning year."""
    demand_fn = Path(
        fn,
        "NT",
        "H2 demand profiles",
        f"H2 {planning_horizon}",
        f"NT_{planning_horizon}.xlsx",
    )
    demand = read_h2_excel(demand_fn, planning_horizon, wscenario)
    demand.columns = [f"{col[:2]} H2" for col in demand.columns]

    return demand


def load_h2_demand(fn: str, planning_horizon: int, wscenario: int) -> pd.DataFrame:
    """
    Load hydrogen demand data for a specific weather scenario and planning year.

    This function retrieves hydrogen demand data from a file, either by loading
    the exact year if available or by performing linear interpolation between
    available years. The data is filtered for a specific weather scenario.

    Parameters
    ----------
    fn : str
        Filepath to the hydrogen demand data file.
    planning_horizon : int
        Planning year for which to retrieve hydrogen demand data.
    wscenario : int
        Weather scenario used to filter the demand data.

    Returns
    -------
    pd.DataFrame
        DataFrame containing hydrogen demand data for the specified planning
        year and weather scenario.
    """

    available_years = get_available_years(fn)
    logger.info(f"Available years: {available_years}, Target year: {planning_horizon}")

    # If target year exists in data, load it directly
    if planning_horizon in available_years:
        logger.info(
            f"Year {planning_horizon} found in available data. Loading directly."
        )
        return load_single_year(fn, planning_horizon, wscenario)

    # Target year not available, do linear interpolation
    return interpolate_demand(
        available_years=available_years,
        planning_horizon=planning_horizon,
        load_single_year_func=load_single_year,
        fn=fn,
        wscenario=wscenario,
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
    planning_horizon = int(snakemake.wildcards.planning_horizons)
    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )
    wscenario = snapshots[0].year
    fn = snakemake.input.h2_demand

    # Check if weather scenario is valid
    wscenario = check_wscenario(wscenario)

    # Load demand with interpolation
    logger.info(
        f"Processing H2 demand for target year: {planning_horizon}, "
        f"weather scenario: {wscenario}"
    )
    demand = load_h2_demand(fn, planning_horizon, wscenario)

    # Reindex demand to fit to snapshots
    demand = align_demand_to_snapshots(demand, snapshots)

    # Export to CSV
    demand.to_csv(snakemake.output.h2_demand, index=True)
