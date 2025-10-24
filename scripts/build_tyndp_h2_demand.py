# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script is used to build TYNDP Scenario Building hydrogen demand to be used in the PyPSA-Eur workflow. The `snapshot` year is used as climatic year (`cyear`). For DE and GA, `cyear` must be one of the following years: 1995, 2008 or 2009. For NT, it must be between 1982 and 2019. If the `snapshot` is not one of these years, then the demand is set to 2009 hydrogen demand (2009 being considered as the most representative of the three years).

Depending on the scenario and hydrogen zone, different years are available in the input data. DE and GA are defined for 2030, 2040 and 2050 for Z2, and 2040, hydrogen zone Z1 is defined in 2050 for DE and GA, in 2040 only for GA. NT scenario is only defined for 2030 and 2040 without a split into two hydrogen zones. All the planning years are read at once. Missing years are linearly interpolated between the available years.
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from _helpers import configure_logging, get_snapshots, set_scenario_config

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
        # rename UK in GB
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


def interpolate_demand(
    available_years: list[int], pyear: int, fn: str, scenario: str, cyear: int
) -> pd.DataFrame:
    """Interpolate demand between available years."""

    # Currently only implemented interpolation and not extrapolation
    lower_years = [y for y in available_years if y < pyear]
    upper_years = [y for y in available_years if y > pyear]

    year_lower = max(lower_years)
    year_upper = min(upper_years)

    logger.info(f"Interpolating {pyear} from {year_lower} and {year_upper}")

    df_lower = load_single_year(fn, scenario, year_lower, cyear)
    df_upper = load_single_year(fn, scenario, year_upper, cyear)

    # Check if data was loaded successfully
    if df_lower.empty and df_upper.empty:
        logger.error("Both years failed to load")
        return pd.DataFrame()
    elif df_lower.empty:
        logger.warning(
            f"Year {year_lower} failed to load. Using zeros for interpolation."
        )
        df_lower = pd.DataFrame(0, index=df_upper.index, columns=df_upper.columns)
    elif df_upper.empty:
        logger.warning(
            f"Year {year_upper} failed to load. Using zeros for interpolation."
        )
        df_upper = pd.DataFrame(0, index=df_lower.index, columns=df_lower.columns)

    df_lower_aligned, df_upper_aligned = df_lower.align(
        df_upper, join="outer", axis=1, fill_value=0
    )

    missing_in_lower = df_upper.columns.difference(df_lower.columns)
    missing_in_upper = df_lower.columns.difference(df_upper.columns)

    if len(missing_in_lower) > 0 or len(missing_in_upper) > 0:
        logger.warning(
            f"Column mismatch between {year_lower} and {year_upper}. "
            f"Missing columns filled with zeros. "
            f"Missing in {year_lower}: {list(missing_in_lower)}, "
            f"Missing in {year_upper}: {list(missing_in_upper)}"
        )

    weight = (pyear - year_lower) / (year_upper - year_lower)
    # Perform linear interpolation
    result = df_lower_aligned * (1 - weight) + df_upper_aligned * weight

    return result


def load_h2_demand(fn: str, scenario: str, pyear: int, cyear: int) -> pd.DataFrame:
    """Load hydrogen demand files into dictionary of dataframes. Filter for specific climatic year and format data."""

    available_years = get_available_years(fn, scenario)
    logger.info(
        f"Scenario {scenario}: Available years: {available_years}, Target year: {pyear}"
    )

    # If target year exists in data, load it directly
    if pyear in available_years:
        logger.info(f"Year {pyear} found in available data. Loading directly.")
        return load_single_year(fn, scenario, pyear, cyear)

    # Target year not available, do linear interpolation
    return interpolate_demand(available_years, pyear, fn, scenario, cyear)


def check_cyear(cyear: int, scenario: str) -> int:
    """Check if the climatic year is valid for the given scenario."""

    valid_years = {
        "NT": np.arange(1983, 2018).tolist(),
        "DE": [1995, 2008, 2009],
        "GA": [1995, 2008, 2009],
    }

    if cyear not in valid_years[scenario]:
        logger.warning(
            "Snapshot year doesn't match available TYNDP data. Falling back to 2009."
        )
        cyear = 2009

    return cyear


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
    cyear = 2009  # snapshots[0].year
    fn = snakemake.input.h2_demand

    # Check if climatic year is valid for scenario
    cyear = check_cyear(cyear, scenario)

    # Load demand with interpolation
    logger.info(
        f"Processing H2 demand for scenario: {scenario}, "
        f"target years: {pyear}, weather year: {cyear}"
    )
    demand = load_h2_demand(fn, scenario, pyear, cyear)

    # Export to CSV
    demand.to_csv(snakemake.output.h2_demand, index=True)
