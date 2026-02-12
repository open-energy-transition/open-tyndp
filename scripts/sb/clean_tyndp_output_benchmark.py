# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script cleans the TYNDP market model output data for benchmarking.

Reads TYNDP market model (MM) MMStandardOutputFile xlsx files from both:
- "Yearly Outputs" sheet for power and electricity data
- "Yearly H2 Outputs" sheet for hydrogen-specific data

Note: Currently, only NT scenario processing is supported.
"""

import logging
import re
from pathlib import Path

import country_converter as coco
import pandas as pd

from scripts._helpers import (
    configure_logging,
    convert_units,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


# Mapping from TYNDP market model (MM) technology names to benchmark carrier names
MM_CARRIER_MAPPING = {
    # Power capacity and generation mappings
    "Nuclear": "nuclear",
    "Lignite old 1": "coal + other fossil",
    "Lignite old 2": "coal + other fossil",
    "Lignite new": "coal + other fossil",
    "Lignite CCS": "coal + other fossil",
    "Lignite biofuel": "biofuels",
    "Hard coal old 1": "coal + other fossil",
    "Hard coal old 2": "coal + other fossil",
    "Hard coal new": "coal + other fossil",
    "Hard coal CCS": "coal + other fossil",
    "Hard Coal biofuel": "biofuels",
    "Gas conventional old 1": "methane",
    "Gas conventional old 2": "methane",
    "Gas CCGT old 1": "methane",
    "Gas CCGT old 2": "methane",
    "Gas CCGT new": "methane",
    "Gas CCGT CCS": "methane",
    "Gas CCGT present 1": "methane",
    "Gas CCGT present 2": "methane",
    "Gas OCGT old": "methane",
    "Gas OCGT new": "methane",
    "Gas biofuel": "biofuels",
    "Light oil": "oil",
    "Heavy oil old 1": "oil",
    "Heavy oil old 2": "oil",
    "Light oil biofuel": "biofuels",
    "Heavy oil biofuel": "biofuels",
    "Oil shale old": "oil",
    "Oil shale new": "oil",
    "Oil shale biofuel": "biofuels",
    "Wind Onshore": "wind onshore",
    "Wind Offshore": "wind offshore",
    "Solar (Photovoltaic)": "solar",
    "Solar (Thermal)": "solar thermal",
    "Run-of-River": "hydro (exc. pump storage)",
    "Reservoir": "hydro (exc. pump storage)",
    "Pondage": "hydro (exc. pump storage)",
    "Pump Storage - Open Loop (turbine)": "hydro and pumped storage",
    "Pump Storage - Closed Loop (turbine)": "hydro and pumped storage",
    "Pump Storage - Open Loop (pump)": "hydro and pumped storage (load)",
    "Pump Storage - Closed Loop (pump)": "hydro and pumped storage (load)",
    "Others renewable": "other res",
    "Others non-renewable": "other non-res",
    "Battery Storage discharge (gen.)": "battery discharge",
    "Battery Storage charge (load)": "battery charge (load)",
    "Hydrogen CCGT": "hydrogen",
    "Hydrogen Fuel Cell": "hydrogen",
    "Demand Side Response Explicit": "demand shedding",
    "Demand Side Response Implicit": "demand shedding",
    # Hydrogen_demand
    "Native Demand (excl. H2 storage charge) [GWhH2]": "exogenous demand",
    # "power generation" - could be calculated from power sheet (H2 Fuel Cell + H2 CCGT) * efficiency
    # "e-fuels" - not available in TYNDP market model NT scenario, included in native demand
    # Hydrogen supply from H2 sheet
    "Electrolyser (gen.)": "p2g",
    "Steam methane reformer": "smr (grey) and smr with ccs (blue)",
    # Note: TYNDP market model doesn't distinguish between grey and blue SMR in the output files
    # "Exchanges with non-modeled nodes" → "imports (renewable & low carbon)" (handled separately)
    # NOT available in TYNDP market model (will not appear in output):
    # - "ammonia imports"
    # - "bi product"
    # - "smr with ccs (blue)" - TYNDP market model only has generic SMR
    # - "undefined for generation" . could be calculated from power sheet (H2 Fuel Cell + H2 CCGT) * efficiency
}


# look up dictionary {name of plot: [sheet_name, output_type]}
LOOKUP_TABLES: dict[str, list[str]] = {
    "power_capacity": ["Yearly Outputs", "Installed Capacities [MW]"],
    "power_generation": ["Yearly Outputs", "Annual generation [GWh]"],
    "hydrogen_demand": [
        "Yearly H2 Outputs",
        "Native Demand (excl. H2 storage charge) [GWhH2]",
    ],
    "hydrogen_supply": ["Yearly H2 Outputs", "Annual generation [GWhH2]"],
    # prices
    "prices_electricity": ["Yearly Outputs", "Marginal Cost Yearly Average [€]"],
    "prices_electricity_excl_shed": [
        "Yearly Outputs",
        "Marginal Cost Yearly Average (excl. 3 000 €/MWh) [€]",
    ],
    "prices_H2": ["Yearly H2 Outputs", "Marginal Cost Yearly Average [€/MWhH2]"],
    "prices_H2_excl_shed": [
        "Yearly H2 Outputs",
        "Marginal Cost Yearly Average (excl. 3 000 €/MWhH2) [€/MWhH2]",
    ],
}

# look up dictionary for crossborder exchanges
CROSS_BORDER_DICT: dict[str, str] = {
    "electricity": "Crossborder exchanges",
    "H2": "Crossborder H2 exchanges",
}


def load_crossborder_sheet(
    sheet_name: str,
    filepath: str | Path,
    skiprows: int = 5,
) -> pd.DataFrame:
    df = pd.read_excel(
        filepath,
        sheet_name=sheet_name,
        skiprows=skiprows,
        usecols=lambda x: x not in [0],  # Skip column indice 0
        index_col=[0],
        nrows=6,
        header=None,
    )

    # rename UK -> GB
    df.iloc[-1, :].replace("UK", "GB")

    # set links names as column
    df = df.set_axis(df.iloc[5], axis=1).drop(df.index[-2:])

    # add buses
    df.loc["bus0", :] = df.columns.str.split("->").str[0]
    df.loc["bus1", :] = df.columns.str.split("->").str[1]

    # rename axis for H2 flows to align with electricity
    df = df.rename(index=lambda x: x.replace("H2", ""))

    # rename axis names
    df.columns.rename("Links", inplace=True)
    df.index.rename("Parameter", inplace=True)

    return df


def load_MM_sheet(
    table_name: str,
    filepath: str | Path,
    eu27: list,
    skiprows: int = 5,
) -> pd.DataFrame:
    """
    Read TYNDP market model sheet from xlsx file.

    Annual sheets have the same structure:
    - Rows 1-4: Metadata (Scenario, Simulator, Date, Status)
    - Row 5: Blank
    - Row 6: Headers (Output type, Output type, then country codes)
    - From row 7: Data rows with [output_type, technology, values...]

    Parameters
    ----------
    filepath : str or Path
        Path to the TYNDP market model xlsx file.
    table_name : str
        Name of the table from LOOKUP_TABLES (e.g., "power_capacity").
    eu27 : list
        List of EU27 country codes.
    skiprows : int, default 5
        Number of metadata rows to skip.

    Returns
    -------
    pd.DataFrame
        DataFrame with multi-index rows [output_type, technology] and country columns.
    """
    sheet_name, output_type = LOOKUP_TABLES[table_name]

    df = pd.read_excel(
        filepath,
        sheet_name=sheet_name,
        skiprows=skiprows,
        header=0,
    )

    # Set multi-index
    level0 = df.iloc[:, 0].ffill()
    level1 = df.iloc[:, 1].fillna(level0)

    # Set as multiindex and keep remaining columns
    df.drop(df.columns[:2], inplace=True, axis=1)
    df.index = pd.MultiIndex.from_arrays([level0, level1])
    df.index.names = ["output_type", "carrier"]

    # Rename and group
    df = df.rename(index=MM_CARRIER_MAPPING, level=1).groupby(level=[0, 1]).sum()
    df = df.loc[output_type]

    # Rename column names to country
    df.rename(
        columns=lambda x: x.replace("UK", "GB").replace("IB", "")[:2], inplace=True
    )
    df_ct = df.T.groupby(df.columns).sum().T

    # Only consider EU27 and sum
    df = df.T.groupby(df.columns).sum().reindex(eu27).sum()

    # Adjust unit
    df = df.rename("value").reset_index()
    df["unit"] = re.search(r"\[(.*)]", output_type).group(1).rstrip("H2")
    final = convert_units(df)
    df_ct["unit"] = re.search(r"\[(.*)]", output_type).group(1).rstrip("H2")
    for col in df_ct.columns[df_ct.columns != "unit"]:
        df_ct = convert_units(df_ct, value_col=col)

    final["table"] = table_name
    return final, df_ct


def set_load_sign(
    MM_data: pd.DataFrame,
    key_word: str = "load",
    tables: list = ["power_generation", "hydrogen_supply"],
) -> pd.DataFrame:
    """
    Set negative sign for load values in market model data.

    Identifies carriers containing the specified keyword in their name
    and negates their values to represent consumption/load.

    Parameters
    ----------
    MM_data : pd.DataFrame
        Market model data with columns 'carrier', 'table', and 'value'.
    key_word : str, default "load"
        Keyword to identify load carriers in the carrier column.
    tables : list, default ["power_generation", "hydrogen_supply"]
        List of table names to apply the sign conversion to.

    Returns
    -------
    pd.DataFrame
        DataFrame with negated values for identified load carriers.
    """
    load_i = MM_data[
        (MM_data.carrier.str.contains(key_word)) & MM_data.table.isin(tables)
    ].index
    MM_data.loc[load_i, "value"] *= -1
    return MM_data


def clean_MM_data_for_benchmarking(MM_data: pd.DataFrame) -> pd.DataFrame:
    """
    Clean market model data for benchmarking analysis.

    Performs the following operations:
    - Removes load and storage discharge entries
    - Excludes pumped storage from power generation table
    - Aggregates hydro capacities (combines regular hydro and pumped storage)
    - Renames carrier categories for consistency with benchmark data

    Parameters
    ----------
    MM_data : pd.DataFrame
        Market model data with columns 'carrier', 'table', and 'value'.
        Expected to contain power capacity and generation data.

    Returns
    -------
    pd.DataFrame
        Cleaned DataFrame with aggregated hydro capacities and
        standardized carrier names.
    """
    # remove load and storages
    MM_data = MM_data[~MM_data.carrier.str.contains("load|discharge")]

    # remove pumped storages from generation
    MM_data = MM_data.query(
        "not(table=='power_generation' and carrier=='hydro and pumped storage')"
    )

    # aggregate hydro capacities
    hydro = (
        MM_data.query(
            "table=='power_capacity' and (carrier == 'hydro (exc. pump storage)' or carrier == 'hydro and pumped storage')"
        )
        .sum(numeric_only=True)
        .value
    )
    MM_data.loc[
        MM_data.query(
            "table=='power_capacity' and carrier=='hydro and pumped storage'"
        ).index,
        "value",
    ] = hydro
    MM_data = MM_data.query(
        "not(table=='power_capacity' and carrier=='hydro (exc. pump storage)')"
    )

    # rename other res to small scale res
    MM_data.loc[
        MM_data.query("table=='power_capacity' and carrier=='other res'").index,
        "carrier",
    ] = "small scale res"
    MM_data.loc[
        MM_data.query("table=='power_capacity' and carrier=='other non-res'").index,
        "carrier",
    ] = "chp and small thermal"

    return MM_data


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_tyndp_output_benchmark",
            planning_horizons="2040",
            scenario="NT",
            configfiles="config/test/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    options = snakemake.params["benchmarking"]
    scenario = snakemake.params["scenario"]
    planning_horizon = int(snakemake.wildcards.planning_horizons)

    # currently only implemented for NT
    if scenario != "NT":
        logger.warning(
            "Processing of TYNDP output files currently only implemented for NT scenario. Exporting empty Data Frame"
        )
        pd.DataFrame().to_csv(snakemake.output.benchmarks)

    # EU27 country codes
    cc = coco.CountryConverter()
    eu27 = cc.EU27as("ISO2").ISO2.tolist()

    # TYNDP market model output file
    tyndp_output_file = snakemake.input.tyndp_output_file

    # Plots for which TYNDP market model output files provide data
    tables_to_process = [
        t for t in LOOKUP_TABLES.keys() if t in options["tables"].keys()
    ]

    logger.info(f"Processing tables: {', '.join(tables_to_process)}")

    benchmarks = {}
    df_ct = {}
    for table in tables_to_process:
        benchmarks[table], df_ct[table] = load_MM_sheet(
            table_name=table, filepath=tyndp_output_file, eu27=eu27, skiprows=5
        )

    MM_data = pd.concat(benchmarks).reset_index(drop=True)
    MM_data_ct = pd.concat(df_ct)

    # set negative sign for loads
    MM_data = set_load_sign(MM_data)

    # clean data for benchmarking
    MM_data = clean_MM_data_for_benchmarking(MM_data)

    MM_data["scenario"] = f"TYNDP {scenario}"
    MM_data["year"] = planning_horizon
    MM_data["source"] = "TYNDP 2024 Market Outputs"

    # load crossborder data
    crossborder = {}
    for key in CROSS_BORDER_DICT.keys():
        crossborder[key] = load_crossborder_sheet(
            CROSS_BORDER_DICT[key], tyndp_output_file
        )
    crossborder_agg = pd.concat(crossborder, axis=1)

    # load prices
    prices_ct = {}
    for table in LOOKUP_TABLES.keys():
        if "price" in table:
            _, prices_ct[table] = load_MM_sheet(
                table_name=table, filepath=tyndp_output_file, eu27=eu27, skiprows=5
            )
    prices = pd.concat(prices_ct)

    # Save data
    MM_data.to_csv(snakemake.output.benchmarks, index=False)
    MM_data_ct.to_csv(snakemake.output.benchmarks_ct)
    crossborder_agg.to_csv(snakemake.output.crossborder)
    prices.to_csv(snakemake.output.prices)
