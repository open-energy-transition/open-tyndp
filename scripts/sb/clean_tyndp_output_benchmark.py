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
import numpy as np
import pandas as pd

from scripts._helpers import (
    align_demand_to_snapshots,
    configure_logging,
    convert_units,
    get_snapshots,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


# Mapping from TYNDP market model (MM) technology names to benchmark carrier names
MM_CARRIER_MAPPING = {
    # Power capacity and generation mappings
    "Nuclear": "nuclear",
    "Lignite old 1": "coal + other fossil (incl. biofuels)",
    "Lignite old 2": "coal + other fossil (incl. biofuels)",
    "Lignite new": "coal + other fossil (incl. biofuels)",
    "Lignite CCS": "coal + other fossil (incl. biofuels)",
    "Lignite biofuel": "coal + other fossil (incl. biofuels)",
    "Hard coal old 1": "coal + other fossil (incl. biofuels)",
    "Hard coal old 2": "coal + other fossil (incl. biofuels)",
    "Hard coal new": "coal + other fossil (incl. biofuels)",
    "Hard coal CCS": "coal + other fossil (incl. biofuels)",
    "Hard Coal biofuel": "coal + other fossil (incl. biofuels)",
    "Gas conventional old 1": "methane (incl. biofuels)",
    "Gas conventional old 2": "methane (incl. biofuels)",
    "Gas CCGT old 1": "methane (incl. biofuels)",
    "Gas CCGT old 2": "methane (incl. biofuels)",
    "Gas CCGT new": "methane (incl. biofuels)",
    "Gas CCGT CCS": "methane (incl. biofuels)",
    "Gas CCGT present 1": "methane (incl. biofuels)",
    "Gas CCGT present 2": "methane (incl. biofuels)",
    "Gas OCGT old": "methane (incl. biofuels)",
    "Gas OCGT new": "methane (incl. biofuels)",
    "Gas biofuel": "methane (incl. biofuels)",
    "Light oil": "oil (incl. biofuels)",
    "Heavy oil old 1": "oil (incl. biofuels)",
    "Heavy oil old 2": "oil (incl. biofuels)",
    "Light oil biofuel": "oil (incl. biofuels)",
    "Heavy oil biofuel": "oil (incl. biofuels)",
    "Oil shale old": "oil (incl. biofuels)",
    "Oil shale new": "oil (incl. biofuels)",
    "Oil shale biofuel": "oil (incl. biofuels)",
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
    "Battery Storage discharge (gen.)": "battery",
    "Battery Storage charge (load)": "battery charge (load)",
    "Hydrogen CCGT": "hydrogen-ccgt",
    "Hydrogen Fuel Cell": "hydrogen-fuel-cell",
    "Demand Side Response Explicit": "demand shedding",
    "Demand Side Response Implicit": "demand shedding",
    # Curtailment/"dump energy"
    "Dump energy [GWh]": "curtailment",
    # Hydrogen_demand
    "Native Demand (excl. H2 storage charge) [GWhH2]": "exogenous demand",
    # "power generation" - could be calculated from power sheet (H2 Fuel Cell + H2 CCGT) * efficiency
    # "e-fuels" - not available in TYNDP market model NT scenario, included in native demand
    # Hydrogen supply from H2 sheet
    "Electrolyser (gen.)": "p2g",
    "Steam methane reformer": "smr (grey) and smr with ccs (blue)",
    # Prices
    "Marginal Cost Yearly Average [€]": "AC",
    "Marginal Cost Yearly Average (excl. 3 000 €/MWh) [€]": "AC",
    "Marginal Cost Yearly Average [€/MWhH2]": "H2",
    "Marginal Cost Yearly Average (excl. 3 000 €/MWhH2) [€/MWhH2]": "H2",
    # Note: TYNDP market model doesn't distinguish between grey and blue SMR in the output files
    # "Exchanges with non-modeled nodes" → "imports (renewable & low carbon)" (handled separately)
    # NOT available in TYNDP market model (will not appear in output):
    # - "ammonia imports"
    # - "bi product"
    # - "smr with ccs (blue)" - TYNDP market model only has generic SMR
    # - "undefined for generation" . could be calculated from power sheet (H2 Fuel Cell + H2 CCGT) * efficiency
}


# look up dictionary {name of plot: [sheet_name, output_type]}
LOOKUP_TABLES: dict[str, list[str| list[str]]] = {
    "power_capacity": ["Yearly Outputs", ["Installed Capacities [MW]"]],
    "power_generation": ["Yearly Outputs", ["Annual generation [GWh]", "Dump energy [GWh]"]],
    "electricity_demand": [
        "Yearly Outputs",
        ["Native Demand (excl. Pump load & Battery charge) [GWh]"],
    ],
    "hydrogen_demand": [
        "Yearly H2 Outputs",
        ["Native Demand (excl. H2 storage charge) [GWhH2]"],
    ],
    "hydrogen_supply": ["Yearly H2 Outputs", ["Annual generation [GWhH2]"]],
    # prices
    "electricity_price": ["Yearly Outputs", ["Marginal Cost Yearly Average [€]"]],
    "electricity_price_excl_shed": [
        "Yearly Outputs",
        ["Marginal Cost Yearly Average (excl. 3 000 €/MWh) [€]"],
    ],
    "hydrogen_price": ["Yearly H2 Outputs", ["Marginal Cost Yearly Average [€/MWhH2]"]],
    "hydrogen_price_excl_shed": [
        "Yearly H2 Outputs",
        ["Marginal Cost Yearly Average (excl. 3 000 €/MWhH2) [€/MWhH2]"],
    ],
}

# look up dictionary for crossborder exchanges
CROSS_BORDER_DICT: dict[str, str] = {
    "electricity": "Crossborder exchanges",
    "H2": "Crossborder H2 exchanges",
}


def normalize_direction(df, cols):
    mask = (df["bus0"] > df["bus1"]) & ~df["bus0"].str.startswith("X")
    assignments = {col: np.where(mask, -df[col], df[col]) for col in cols}

    # Add bus0, bus1, and border to assignments
    assignments.update(
        {
            "bus0": np.where(mask, df["bus1"], df["bus0"]),
            "bus1": np.where(mask, df["bus0"], df["bus1"]),
            "border": lambda df: df.bus0 + "->" + df.bus1,
        }
    )

    return df.assign(**assignments)


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

    # set links names as column
    df = df.set_axis(df.iloc[5], axis=1).drop(df.index[-2:])

    # Rename column names
    df.rename(columns=lambda x: x.replace("UK", "GB"), inplace=True)
    df.rename(columns=lambda x: x.replace("_", " "), inplace=True)  # for H2

    # normalize direction
    attributes = df.index
    # add buses
    df.loc["bus0", :] = df.columns.str.split("->").str[0]
    df.loc["bus1", :] = df.columns.str.split("->").str[1]
    df = normalize_direction(df.T.reset_index(drop=True), cols=attributes)

    # set index
    df = df.set_index("border").T
    df.index.rename("Parameter", inplace=True)
    # rename axis for H2 flows to align with electricity
    df = df.rename(index=lambda x: x.replace("H2", ""))

    # convert units
    df.loc["Sum [MWh]:"] = df.loc["Sum [GWh]:"].astype(float) * 1e3
    df.drop("Sum [GWh]:", inplace=True)
    df["unit"] = df.index.str.extract(r"\[(.*?)\]", expand=False)

    # rename index
    stats_labels = {
        "Avg [MW]:": "avg",
        "Max [MW]:": "max",
        "Min [MW]:": "min",
        "Sum [MWh]:": "sum",
    }
    df.rename(index=stats_labels, inplace=True)

    return df.sort_index()


def load_MM_sheet(
    table_name: str,
    filepath: str | Path,
    countries: list[str],
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
    countries : list[str]
        List of modelled countries
    eu27 : list
        List of EU27 country codes.
    skiprows : int, default 5
        Number of metadata rows to skip.

    Returns
    -------
    pd.DataFrame
        Market Model data in long format (incl. EU27).
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
    df = df.loc[output_type].droplevel(0)

    # Rename and filter column names (buses)
    df.rename(
        columns=lambda x: x.replace("UK", "GB").replace("_H2", " H2"),
        inplace=True,
    )
    op = "sum" if "price" not in table_name else "mean"
    df_nodal = (
        df.T.groupby(df.columns)
        .agg(op)
        .T.reset_index()
        .melt(id_vars=["carrier"], var_name="bus")
    )
    df_nodal = df_nodal[df_nodal.bus.str.extract(r"^(?:IB)?(.{2})")[0].isin(countries)]

    # Add EU27 (load-weighted average for prices)
    df_eu27 = df_nodal[df_nodal.bus.str.extract(r"^(?:IB)?(.{2})")[0].isin(eu27)]

    if "price" in table_name:
        weights = (
            load_MM_sheet(
                table_name=f"{table_name.split('_')[0]}_demand",
                filepath=tyndp_output_file,
                countries=countries,
                eu27=eu27,
                skiprows=5,
            )
            .set_index("bus")
            .value
        )
        normalizer = weights.sum()
    else:
        weights = pd.Series(1.0, index=df_eu27.bus.unique())
        normalizer = 1

    df_eu27 = (
        df_eu27.assign(value=lambda x: x.bus.map(weights).fillna(0) * x.value)
        .groupby(by=["carrier"])
        .value.sum()
        .div(normalizer)
        .reset_index()
        .assign(bus="EU27")
    )
    df = pd.concat([df_nodal, df_eu27])

    # Add metadata
    df["unit"] = re.sub(
        r"€(/MWh)?",
        "EUR/MWh",
        re.search(r"\[(.*)]", output_type[0]).group(1).rstrip("H2"),
    )

    df["table"] = table_name
    if "price" not in table_name:
        df = convert_units(df)

    return df


def load_h2_demand_ts(
    sheet_name: str, filepath: str | Path, snapshots: pd.DatetimeIndex
) -> pd.DataFrame:
    """
    Load hourly H2 demand time series from a TYNDP Market Model Outputs Excel file.

    Parameters
    ----------
    sheet_name : str
        Name of the Excel sheet containing hourly H2 data.
    filepath : str or Path
        Path to the TYNDP Market Model Outputs Excel file.
    snapshots : pd.DatetimeIndex
        Model snapshot index.

    Returns
    -------
    pd.DataFrame
        DataFrame of hourly H2 demand.
    """
    df = (
        pd.read_excel(
            filepath,
            sheet_name=sheet_name,
            skiprows=12,
            index_col=1,
        )
        .filter(like="H2_LOAD")
        .rename(lambda s: s.split("_")[0].replace("UK", "GB") + " H2", axis=1)
    )

    df = align_demand_to_snapshots(df, snapshots, format="%d%b%H:%M")

    return df


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


def clean_MM_data_for_benchmarking(
    MM_data: pd.DataFrame, offshore_hubs: bool = False
) -> pd.DataFrame:
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
    offshore_hubs : bool, default False
        Whether offshore hubs are modeled.

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

    # remove batteries from generation
    MM_data = MM_data.query("not(table=='power_generation' and carrier=='battery')")

    # aggregate hydro capacities
    mask = (MM_data.table == "power_capacity") & (
        MM_data.carrier == "hydro (exc. pump storage)"
    )
    MM_data.loc[mask, "carrier"] = "hydro and pumped storage"
    MM_data = (
        MM_data.groupby([c for c in MM_data.columns if c != "value"])
        .sum()
        .reset_index()
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

    # reflect H2 CCGT and fuel cells consumptions in the yearly H2 demand
    h2_power = MM_data.query(
        "table=='power_generation' and carrier.isin(['hydrogen-ccgt', 'hydrogen-fuel-cell'])"
    ).copy()
    eff_fuel_cell, eff_ccgt = 0.5, 0.59  # TODO Remove hard coded values
    mask_eu27 = ~(h2_power["bus"] == "EU27")
    h2_power.loc[mask_eu27, "bus"] = h2_power.loc[mask_eu27, "bus"].str[:2] + " H2"
    h2_power.loc[h2_power.carrier == "hydrogen-ccgt", "value"] = h2_power.loc[
        h2_power.carrier == "hydrogen-ccgt", "value"
    ].div(eff_ccgt)
    h2_power.loc[h2_power.carrier == "hydrogen-fuel-cell", "value"] = h2_power.loc[
        h2_power.carrier == "hydrogen-fuel-cell", "value"
    ].div(eff_fuel_cell)
    h2_power.loc[:, "table"] = "hydrogen_demand"
    h2_power.loc[:, "carrier"] = "power generation"

    MM_data = (
        pd.concat([MM_data, h2_power])
        .replace({"hydrogen-ccgt": "hydrogen", "hydrogen-fuel-cell": "hydrogen"})
        .groupby([c for c in MM_data.columns if c != "value"])
        .sum()
        .reset_index()
    )

    # Norway (NO) has no H2 demand, so its reported H2 price is set to 0 EUR/MWh_H2
    # to avoid misleading benchmark errors
    MM_data.loc[
        MM_data.table.str.contains("hydrogen_price") & MM_data.bus.str.contains("NO"),
        "value",
    ] = 0

    # When offshore hubs are modeled, exclude price data for conventional offshore nodes
    if offshore_hubs:
        offshore_node_conv = [  # noqa: F841
            "BEOF",
            "DEKF",
            "DKBH",
            "DKKF",
            "DKNS",
            "EEOF",
            "LTOF",
            "NL60",
            "NL6H",
            "NLA0",
            "NLBH",
            "NLLL",
        ]
        MM_data = MM_data.query(
            "~(bus.isin(@offshore_node_conv) and table.str.contains('price'))"
        )

    return MM_data


def assign_meta_data(df, planning_horizon, scenario):
    df["scenario"] = f"TYNDP {scenario}"
    df["year"] = planning_horizon
    df["source"] = "TYNDP 2024 Market Model Outputs"


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
    countries = snakemake.params["countries"]

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
    for table in tables_to_process:
        benchmarks[table] = load_MM_sheet(
            table_name=table,
            filepath=tyndp_output_file,
            countries=countries,
            eu27=eu27,
            skiprows=5,
        )

    MM_data = pd.concat(benchmarks).reset_index(drop=True)

    # set negative sign for loads
    MM_data = set_load_sign(MM_data)

    # clean data for benchmarking
    MM_data = clean_MM_data_for_benchmarking(
        MM_data, offshore_hubs=snakemake.params.offshore_hubs
    )

    # load crossborder data
    logger.info(
        f"Processing tables of cross-border flows for: {', '.join(CROSS_BORDER_DICT.keys())}"
    )
    crossborder = {}
    for key in CROSS_BORDER_DICT.keys():
        crossborder[key] = load_crossborder_sheet(
            CROSS_BORDER_DICT[key], tyndp_output_file
        )
    crossborder_agg = pd.concat(crossborder, axis=1)

    # load h2 demand time series
    logger.info("Processing hourly H2 demand tables")
    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )
    h2_demand_ts = load_h2_demand_ts(
        sheet_name="Hourly H2 Data", filepath=tyndp_output_file, snapshots=snapshots
    )

    # assign meta data
    assign_meta_data(MM_data, planning_horizon, scenario)
    assign_meta_data(crossborder_agg, planning_horizon, scenario)
    assign_meta_data(h2_demand_ts, planning_horizon, scenario)

    # Save data
    MM_data.to_csv(snakemake.output.benchmarks, index=False)
    crossborder_agg.to_csv(snakemake.output.crossborder)
    h2_demand_ts.to_csv(snakemake.output.h2_demand)
