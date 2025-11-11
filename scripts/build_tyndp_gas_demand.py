# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP Scenario Building gas demand for Open-TYNDP.

This script processes methane (gas) demand data from TYNDP 2024 Supply Tool,
extracting both final energy demand (direct gas consumption) and heat-related
gas demand. The data is filtered and interpolated based on the selected
scenario (Distributed Energy, Global Ambition, or National Trends) and
planning horizon.

Data Availability by Scenario
------------------------------

The input data has different temporal coverage depending on scenario:

**National Trends (NT)**:
  - Available for 2030 and 2040
  - Includes final gas demand (incl. heat demand) from NT+ data collection
  - Heat demand processed with distribution shares and efficiency factors
  - Special handling for Italian demand

**Distributed Energy (DE) and Global Ambition (GA)**:
  - Processing not yet implemented

Processing
----------

Missing years are linearly interpolated between available data points.

For each planning year, the script:
1. Reads final energy demand (direct gas consumption by carrier)
2. Reads heat demand and converts it to primary energy using distribution and efficiency data
3. Combines both demand components
4. Aggregates by country/bus

Inputs
------

- ``data/tyndp_2024_bundle/Supply Tool/20240518-Supply-Tool.xlsm``: TYNDP 2024 Supply Tool Excel file containing:
  - NT+ data sheet: Final demand by country
  - Other data and Conversions sheet: Heat distribution and efficiency factors
  - IT sheet: Italian gas production data

Outputs
-------

- ``gas_demand_tyndp_{planning_horizons}.csv``: Processed gas demand data (MWh) by country/bus
  for the specified planning horizon
"""

import logging
from typing import Literal

import country_converter as coco
import pandas as pd
from _helpers import configure_logging, interpolate_demand, set_scenario_config

logger = logging.getLogger(__name__)
cc = coco.CountryConverter()


def convert_country_to_bus(names: pd.Index) -> pd.Index:
    """Convert country names to bus names."""
    names = pd.Index(cc.convert(names, to="iso2"))
    return names + "00"


def read_fed_data(fn: str, scenario: str, pyear: int) -> tuple[pd.Series, pd.Series]:
    """
    Read and process final gas demand data and final heat demand data from Supply Tool for a specific year.
    """
    try:
        demand_fed = pd.read_excel(
            fn,
            usecols="B:AC",
            header=1,
            index_col=0,
            nrows=25,
            skiprows=27 if pyear == 2040 else 0,
            sheet_name="NT+ data",
        )

        # Set buses as column names
        demand_fed.columns = convert_country_to_bus(demand_fed.columns)

        # Extract final heat demand
        demand_heat = demand_fed.loc["Heat"].mul(1e3)  # MWh

        # Extract final methane demand
        demand_fed = (
            demand_fed.loc[
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

    except Exception as e:
        logger.warning(
            f"Failed to read final gas demand for scenario {scenario} and pyear {pyear}: "
            f"{type(e).__name__}: {e}"
        )
        demand_fed = pd.Series()
        demand_heat = pd.Series()

    return demand_fed, demand_heat


def read_heat_frame(
    fn: str, pyear: int, type: Literal["distribution", "efficiency"]
) -> pd.DataFrame:
    """Read heat distribution and efficiency tables in 'Other data and Conversions' sheet of Supply Tool."""
    if type not in ["distribution", "efficiency"]:
        raise ValueError(
            f"Invalid type '{type}'. Must be 'distribution' or 'efficiency'."
        )
    if pyear not in [2030, 2040, 2050]:
        raise ValueError(f"Invalid pyear '{pyear}'. Must be 2030, 2040, or 2050.")

    offset = ((pyear - 2030) // 10) * 34
    offset += 17 if type == "efficiency" else 0

    df = pd.read_excel(
        fn,
        usecols="A:AB",
        header=1,
        index_col=0,
        nrows=16,
        skiprows=45 + offset,
        sheet_name="Other data and Conversions ",
    )

    # Set buses as column names
    df.columns = convert_country_to_bus(df.columns)

    return df


def read_it_gas_prod(fn: str, pyear: int) -> float:
    """Read Italian gas production used for heat. This data is hardcoded in Supply Tool IT sheet."""
    return (
        pd.read_excel(
            fn,
            usecols="H:K",
            header=0,
            index_col=0,
            nrows=2,
            skiprows=31,
            sheet_name="IT",
        ).loc[pyear, "For heat"]
        * 1e6
    )  # MWh


def read_heat_data(
    heat_fed: pd.Series, fn: str, scenario: str, pyear: int
) -> pd.Series:
    """Read and process heat-related gas demand data from Supply Tool for a specific year."""
    try:
        shares = read_heat_frame(fn, pyear, "distribution")
        efficiencies = read_heat_frame(fn, pyear, "efficiency")

        # Compute heat primary energy demand
        demand_primary = heat_fed * shares * (1 / efficiencies)

        # Filter energy carriers
        demand = demand_primary.loc[
            [
                "Biogas",
                "E-Methane",
                "Natural gas",
                "Other fossil gas",
                "Waste gas",
            ]
        ].sum()  # MWh

        # Apply Italian solution as specified in Supply Tool
        demand.loc["IT00"] = read_it_gas_prod(fn, pyear)

    except Exception as e:
        logger.warning(
            f"Failed to read heat demand data for scenario {scenario} and pyear {pyear}: "
            f"{type(e).__name__}: {e}"
        )
        demand = pd.Series()

    return demand


def read_supply_tool(fn: str, scenario: str, pyear: int) -> pd.Series:
    """Read and process both final gas demand and heat demand data from Supply Tool for a specific year."""
    demand_fed, heat_fed = read_fed_data(fn, scenario, pyear)
    demand_heat = read_heat_data(heat_fed, fn, scenario, pyear)

    demand = pd.concat([demand_fed, demand_heat], axis=1).sum(axis=1)
    demand.name = "p_nom"

    return demand


def load_single_year(fn: str, scenario: str, pyear: int) -> pd.Series:
    """Load demand data for a single planning year."""
    if scenario == "NT":
        demand = read_supply_tool(fn, scenario, pyear)
    elif scenario in ["DE", "GA"]:
        # TODO Implement processing for DE/GA
        demand = pd.Series()

    return demand


def load_gas_demand(fn: str, scenario: str, pyear: int) -> pd.Series:
    """
    Load gas demand data for a specific scenario and planning year.

    This function retrieves gas demand data from a file, either by loading
    the exact year if available or by performing linear interpolation between
    available years.

    Parameters
    ----------
    fn : str
        Filepath to the gas demand data file.
    scenario : str
        Name of the scenario to load.
    pyear : int
        Planning year for which to retrieve gas demand data.

    Returns
    -------
    pd.Series
        Series containing gas demand data for the specified scenario and planning year.
    """

    available_years = [2030, 2040]

    # If target year exists in data, load it directly
    if pyear in available_years:
        logger.debug(f"Year {pyear} found in available data. Loading directly.")
        return load_single_year(fn, scenario, pyear)

    # Target year not available, do linear interpolation
    return interpolate_demand(
        available_years=available_years,
        pyear=pyear,
        load_single_year_fun=load_single_year,
        fn=fn,
        scenario=scenario,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_gas_demand",
            configfiles="config/test/config.tyndp.yaml",
            planning_horizons=2035,
            run="NT",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    scenario = snakemake.params["scenario"]
    fn = snakemake.input.supply_tool
    pyear = int(snakemake.wildcards.planning_horizons)

    if scenario != "NT":
        logger.warning(f"Gas demand processing is not supported yet for {scenario}.")

    # Load demand with interpolation
    logger.info(f"Processing gas demand for scenario: {scenario}")
    demand = load_gas_demand(fn, scenario, pyear)

    # Export to CSV
    demand.to_csv(snakemake.output.gas_demand, index=True)
