# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP Scenario Building gas demand for Open-TYNDP.

# TODO Update the documentation

This script processes methane demand data from TYNDP 2024, using the
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
from bisect import bisect_right

import country_converter as coco
import pandas as pd
from _helpers import (
    configure_logging,
    set_scenario_config,
)

logger = logging.getLogger(__name__)
cc = coco.CountryConverter()


def read_supply_tool(fn: str, scenario: str, pyear: int) -> pd.DataFrame:
    """Read and process gas demand data from Supply Tool for a specific year."""
    try:
        demand = pd.read_excel(
            fn,
            usecols="B:AC",
            header=1,
            index_col=0,
            nrows=25,
            skiprows=27 if pyear == 2040 else 0,
            sheet_name="NT+ data",
        )

        # Set nodes as column name
        demand.columns = cc.convert(demand.columns, to="iso2")
        demand.columns = demand.columns + "00"

        # Filter energy carriers
        demand = (
            demand.loc[
                [
                    "E-Methane",
                    "Other fossil gas",
                    "Biomethane",
                    "Natural gas",
                    "Waste gas",
                    "Gas for Cooking",
                    "Methane (LNG)",
                ]
            ]
            .mul(1e3)
            .sum()
        )  # MWh
        demand.name = "p_nom"

    except Exception as e:
        logger.warning(
            f"Failed to read gas demand for scenario {scenario} and pyear {pyear}: "
            f"{type(e).__name__}: {e}"
        )
        demand = pd.DataFrame()

    return demand


def load_single_year(fn: str, scenario: str, pyear: int) -> pd.DataFrame:
    """Load demand data for a single planning year."""
    if scenario == "NT":
        demand = read_supply_tool(fn, scenario, pyear)
    elif scenario in ["DE", "GA"]:
        # ToDo Implement processing for DE/GA
        pass

    return demand


def interpolate_demand(
    available_years: list[int], pyear: int, fn: str, scenario: str
) -> pd.DataFrame:
    """Interpolate demand between available years."""

    # Currently only implemented interpolation and not extrapolation
    idx = bisect_right(available_years, pyear)
    if idx == 0:
        # Planning horizon is before all available years
        logger.warning(
            f"Year {pyear} is before the first available year {available_years[0]}. "
            f"Falling back to first available year."
        )
        year_lower = year_upper = available_years[0]
    elif idx == len(available_years):
        # Planning horizon is after all available years
        logger.warning(
            f"Year {pyear} is after the latest available year {available_years[-1]}. "
            f"Falling back to latest available year."
        )
        year_lower = year_upper = available_years[-1]
    else:
        year_lower = available_years[idx - 1]
        year_upper = available_years[idx]

    logger.info(f"Interpolating {pyear} from {year_lower} and {year_upper}")

    df_lower = load_single_year(fn, scenario, year_lower)
    df_upper = load_single_year(fn, scenario, year_upper)

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
            f"Year {year_upper} failed to load. Using data from lower year for interpolation."
        )
        df_upper = df_lower

    missing_in_lower = df_upper.columns.difference(df_lower.columns)
    missing_in_upper = df_lower.columns.difference(df_upper.columns)

    if len(missing_in_lower) > 0 or len(missing_in_upper) > 0:
        logger.warning(
            f"Column mismatch between {year_lower} and {year_upper}. "
            f"Missing columns filled with zeros. "
            f"Missing in {year_lower}: {list(missing_in_lower)}, "
            f"Missing in {year_upper}: {list(missing_in_upper)}"
        )
    df_lower_aligned, df_upper_aligned = df_lower.align(
        df_upper, join="outer", axis=1, fill_value=0
    )

    weight = (pyear - year_lower) / (year_upper - year_lower)
    # Perform linear interpolation
    result = df_lower_aligned * (1 - weight) + df_upper_aligned * weight

    return result


def load_gas_demand(fn: str, scenario: str, pyear: int) -> pd.DataFrame:
    """
    Load gas demand data for a specific scenario and planning year.

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

    Returns
    -------
    pd.DataFrame
        DataFrame containing gas demand data for the specified scenario and planning year.
    """

    available_years = [2030, 2040]

    # If target year exists in data, load it directly
    if pyear in available_years:
        logger.info(f"Year {pyear} found in available data. Loading directly.")
        return load_single_year(fn, scenario, pyear)

    # Target year not available, do linear interpolation
    return interpolate_demand(available_years, pyear, fn, scenario)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_gas_demand",
            configfiles="config/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    scenario = snakemake.params["scenario"]
    fn = snakemake.input.supply_tool
    pyear = int(snakemake.wildcards.planning_horizons)

    if scenario != "NT":
        logger.warning(
            f"Gas demand processing is not supported yet for {scenario}. Falling back to NT data."
        )

    # Load demand with interpolation
    logger.info(f"Processing gas demand for scenario: {scenario}")
    demand = load_gas_demand(fn, scenario, pyear)

    # Export to CSV
    demand.to_csv(snakemake.output.gas_demand, index=True)
