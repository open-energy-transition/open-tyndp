# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script cleans the TYNDP market model output data for benchmarking.

Reads TYNDP market model (MM) TimeSeries Dashboard xlsx files, one per country, 
per planning horizon and weather scenario
- "Installed Capacity" sheet for installed capacities per market zone
- Nodal sheets for aggregates of hourly generation, load, prices, curtailment and unserved energy
- "Exchanges" sheet for cross-border electricity and hydrogen flows

Note: Currently, only NT scenario processing is supported.
"""

import logging
import re
from pathlib import Path

import country_converter as coco
import numpy as np
import pandas as pd
from functools import partial

from scripts._helpers import (
    align_demand_to_snapshots,
    configure_logging,
    convert_units,
    get_snapshots,
    normalize_direction,
    set_scenario_config,
    format_bz_names,
    get_weather_scenario,
)

from scripts.build_tyndp_network import extract_country


logger = logging.getLogger(__name__)

# List of sheets for both electricity and hydrogen timeseries
ELECTRICITY_SHEETS = ["E-Market", "Prosumer", "Offshore"]
H2_SHEETS = ["H2 Zone 1", "H2 Zone 2"]

# look up dictionary {name of plot: [sheet_name, output_type]}
LOOKUP_TABLES: dict[str, dict] = {
    "power_capacity": {"sheet": ["Installed Capacity"]},
    "power_generation": {"sheet": ELECTRICITY_SHEETS, "stats": "sum"},
    "electricity_demand": {
        "sheet": ELECTRICITY_SHEETS, 
        "category":["Native Demand [MW_e]", "Fixed Demand [MW_e]"], 
        "stats": "sum", 
    },
    "hydrogen_demand": {
        "sheet": H2_SHEETS, "category":["Native Demand [MW_H2]"], "stats": "sum", 
    },
    "hydrogen_supply": {"sheet": H2_SHEETS, "stats": "sum"},
    # TODO: For shedding hours, TYNDP 2026 does not give a direct value, so it can be derived
    # where "Energy Not Served" is above 0.
    # "electricity_demand_shedding_hours": [
        #"Yearly Outputs",
        #"Loss of load expectation [hour]  ",
    #],  # includes white space
    #"hydrogen_demand_shedding_hours": [
        #"Yearly H2 Outputs",
        #"Loss of H2 load expectation [hour]  ",
    #],  # includes white space
    # prices
    "electricity_price": {
        "sheet": ["E-Market"], "category":["Marginal Cost [€/MWh_e]"], "stats": "avg", 
    },
    #"electricity_price_excl_shed": [
       # "Yearly Outputs",
        #"Marginal Cost Yearly Average (excl. 3 000 €/MWh) [€]",
    #],
    "hydrogen_price": {
        "sheet": ["H2 Zone 2"], "category":["Marginal Cost [€/MWh_H2]"], "stats": "avg", 
    },
    #"hydrogen_price_excl_shed": [
        #"Yearly H2 Outputs",
        #"Marginal Cost Yearly Average (excl. 3 000 €/MWhH2) [€/MWhH2]",
    #],
}

# look up dictionary for crossborder exchanges
CROSS_BORDER_DICT: dict[str, str] = {
    "electricity": "E-Market Exchanges [MW_e]",
    "H2": "H2 Zone 2 Exchanges [MW_H2]",
    "H2 imports": "H2 Imports to H2 Zone 2 [MW_H2]",
}


def _load_mm_carrier_mapping(
    carrier_mapping_fn: str, tables: dict
) -> tuple[dict[str, dict], dict[str, dict]]:
    """
    Load mapping from TYNDP Market Model (MM) carrier names to benchmark carrier names.
    """
    tech_map = pd.read_csv(carrier_mapping_fn)
    output_map = {}
    for table, table_opts in tables.items():
        if table not in LOOKUP_TABLES:
            continue
        col = table_opts.get("mapping_col")
        if col is None:
            continue
        if col not in tech_map.columns:
            logger.warning(
                f"No existing mapping for table {table} in 'tyndp_technology_map.csv'."
            )
            continue
        output_map[table] = (
            tech_map[["tyndp_output_carrier", col]]
            .dropna(subset=["tyndp_output_carrier", col])
            .set_index("tyndp_output_carrier")[col]
            .to_dict()
        )

    return output_map

def split_unit_from_category(carrier_labels: pd.Series) -> tuple[pd.Series, pd.Series]:
    """
    Split units from category (carrier) labels in the dashboard
    """
    unit = (
        carrier_labels.str.extract(r"\[([^\]]*)\]", expand=False)
        .str.replace("€", "EUR", regex=False)
        .str.replace(r"_(e|H2)$", "", regex=True)
    )
    return carrier_labels.str.replace(r"\s*\[[^\]]*\]\s*$", "", regex=True).str.strip(), unit

def rename_prosumer_nodes(buses:pd.Series) -> pd.Series:
    """
    Rename prosumer nodes to match Open-TYNDP bus names 
    """
    return format_bz_names(buses).str.removesuffix("RETE")

def load_dashboard_sheet(
    filepath: str | Path,
    sheet_name: str,
    stats: list[str],
    categories: list[str] = None,
) -> pd.DataFrame:
    """
    Load a sheet from a TYNDP TimeSeries Dashboard output file.

    All sheets, except "Installed Capacity", share the same structure:
    - Rows 1-4: Yearly values (min, max, average, and sum) of the hourly timeseries
    - Row 6: Market Zone
    - Row 7: Category (similar to carrier)
    - Row 8: Element (similar to type)
    - Row 9: Generation/Load, indicates flow direction.

    Parameters
    ----------
    filepath : str or Path
        Path to the Excel file
    sheet_name : str
        Name of the Excel sheet to read
    stats : list[str]
        Yearly stats to extract, a subset of ["min", "max", "avg", "sum"]

    Returns
    -------
    pd.DataFrame
        Long format with columns [sheet, bus, carrier, element, flow, stats,
        unit, value]
    """
    stats_rows = {"min": 0, "max": 1, "avg": 2, "sum": 3}
    row_zone, row_carrier, row_element, row_flow, col_data = 5, 6, 7, 8, 2

    # Load the file and check if the sheet exists as some nodes are missing sheets.
    file = pd.ExcelFile(filepath, engine="calamine")
    if sheet_name not in file.sheet_names:
        return pd.DataFrame()

    df = pd.read_excel(
        file, sheet_name=sheet_name, header=None, nrows=row_flow + 1
    )
    cols = df.columns[col_data:]
    if categories: 
        cols = cols[df.loc[row_category cols].isin(categories)]

    carrier, unit = split_unit_from_category(df.iloc[row_category,cols:].astype(str))
    df = pd.Dataframe(
        {
            "sheet": sheet_name,
            "bus": df.loc[row_zone, cols],
            "carrier": carrier,
            "element": df.loc[row_element, cols],
            "flow": df.loc[row_flow, cols],
            "unit": unit,
        }
    )

    return pd.concat()


def load_installed_capacity(filepath: str | Path,) -> pd.DataFrame:
    """
    Load the sheet "Installed Capacity" from a TYNDP TimeSeries Dashboard output file.

    Parameters
    ----------
    filepath : str or Path
        Path to the Excel file

    Returns
    -------
    pd.DataFrame
        Long format with columns [sheet, bus, carrier, unit, value]
    """
    row_zone, row_data = 1, 2

    df = pd.read_excel(
        filepath, sheet_name="Installed Capacity", header=None, engine="calamine"
    )
    carrier, unit = split_unit_from_category(df.iloc[row_data:, 0].astype(str))

    return( 
        df.iloc[row_data:, 1:]
        .apply(pd.to_numeric, errors="coerce")
        .set_axis(pd.Index(df.iloc[row_zone, 1:], name="bus"), axis=1)
        .set_axis(pd.MultiIndex.from_arrays([carrier, unit], names=["carrier", "unit"]))
        .stack(future_stack=True)
        .rename("value")
        .reset_index()
        .dropna(subset=["value"])
        .assign(sheet="Installed Capacity")
    )


def load_crossborder_sheet(
    sheet_name: str,
    filepath: list[Path],
    stats: list[str] = ["min", "max", "avg", "sum"],
) -> pd.DataFrame:
    """
    Load the cross-border flow sheet from a TYNDP Market Model output file.

    In TYNDP 2026, both electricity and hydrogen carriers share the "Exchanges" sheet.

    Parameters
    ----------

    sheet_name : str
        Name of the Excel sheet to read
    filepath : list[Path]
        Paths to the Excel files
    stats : list[str], optional
        Yearly aggregates to extract

    Returns
    -------
    pd.DataFrame
        DataFrame with normalized cross-border flow data.
    """
    df = pd.concat(
        [load_dashboard_sheet(sheet_name,fn, stats) for fn in filepath]
    )

    categories = {v: k for k, v in CROSS_BORDER_DICT.items()}

    df = (
        df.assign(
            carrier=lambda x: x.carrier.map(categories),
            border=lambda x: format_bz_names(x.element),
        )
        .drop(subset=["carrier"])
        .drop_duplicates(["border", "stats"])
        .pivot(index=["carrier", "bus", "border"], columns="stats", values="value")
        .reset_index("carrier", "bus")
    )

    # normalize direction
    df = normalize_direction(df, cols=stats, buses_from_index=True, connector="-")
    mask = df["min"] > df["max"]
    df.loc[mask, ["min", "max"]] = df.loc[mask, ["max", "min"]].values

    # convert units
    df["sum"] = df["sum"].mul(1e3)

    return df.sort_index()


def load_MM_sheet(
    table_name: str,
    filepath: str | Path,
    countries: list[str],
    eu27: list,
    mapping: dict[str, dict[str, str]],
) -> pd.DataFrame:
    """
    Read benchmarking table from TYNDP 2026 market model output files
    

    Parameters
    ----------
    filepath : str or Path
        Path to the TYNDP market model xlsx file.
    table_name : str
        Name of the table from LOOKUP_TABLES (e.g., "power_capacity").
    countries : list[str]
        List of modelled countries.
    eu27 : list
        List of EU27 country codes.
    mapping : dict[str, dict[str, str]]
        Carrier mapping from market model carrier names to benchmarking carrier names per table.

    Returns
    -------
    pd.DataFrame
        Market Model data in long format (incl. EU27).
    """
    opt = LOOKUP_TABLES[table_name]

    if table_name == "power_capacity":
        df = pd.concat([load_installed_capacity(fn) for fn in filepath])
    else:
        df = pd.concat(
            [
                load_dashboard_sheet(fn, sheet, [opt["stat"]], opt.get("categories)"))
                for sheet in opt["sheets"]
                for fn in filepath
            ]
        )

    # Only include mapped carriers
    carriers = df.carrier.map(mapping[table_name])
    if unmapped := sorted(df.carrier[carriers.isna()].unique()):
        logger.warning(
            f"No carrier mappings for table '{table_name}' for: {unmapped}. They will be excluded from the benchmarking for this table."
        )
    df = df.assign(carrier=carriers).dropna(subset=["carrier"])

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

    # Add EU27 / Pan-EU load-weighted average for prices
    if "price" in table_name:
        df_eu = df_nodal
        bus_name = "Pan-EU"
        weights = (
            load_MM_sheet(
                table_name=f"{table_name.split('_')[0]}_demand",
                filepath=tyndp_output_file,
                countries=countries,
                eu27=eu27,
                mapping=mapping,
                skiprows=5,
            )
            .query("bus!='EU27'")
            .set_index("bus")
            .value
        )
        normalizer = weights.sum()
    else:
        df_eu = df_nodal[df_nodal.bus.str.extract(r"^(?:IB)?(.{2})")[0].isin(eu27)]
        bus_name = "EU27"
        weights = pd.Series(1.0, index=df_eu.bus.unique())
        normalizer = 1

    df_eu = (
        df_eu.assign(value=lambda x: x.bus.map(weights) * x.value)
        .groupby(by=["carrier"])
        .value.sum()
        .div(normalizer)
        .reset_index()
        .assign(bus=bus_name)
    )
    df = pd.concat([df_nodal, df_eu])

    df["table"] = table_name
    if "price" not in table_name and "hours" not in table_name:
        df = convert_units(df)

    return df


def load_demand_ts(
    sheet_name: list[str],
    filepath: list[Path],
    snapshots: pd.DatetimeIndex,
    categories: list[str],
) -> pd.DataFrame:
    """
    Load hourly demand time series from a TYNDP Market Model Outputs Excel file.

    Parameters
    ----------
    sheet_name : str
        Name of the Excel sheet containing hourly data.
    filepath : str or Path
        Path to the TYNDP Market Model Outputs Excel file.
    snapshots : pd.DatetimeIndex
        Model snapshot index.
    carrier : str
        Energy carrier to load.

    Returns
    -------
    pd.DataFrame
        DataFrame of hourly demand indexed by snapshot.
    """
    row_zone, row_category, row_data, col_time = 5, 6, 10, 1

    dfs =[]
    for filepath in filepath:
        file = pd.ExcelFile(filepath, engine="calamine")
        for sheet_name in [s for s in sheet_name if s in file.sheet_name]:
            df = pd.read_excel(file, sheet_name=sheet_name, header=None)
            cols = df.columns[2:]
            cols = cols[df.loc[row_category, cols].isin(categories)]

            demand = df.loc[row_data:, cols].apply(pd.to_numeric, errors="coerce")
            demand.columns = rename_prosumer_nodes(df.loc[row_zone, cols].astype(str))
            demand.index = df.loc[row_data:, col_time]
            dfs.append(demand)

    # Prosumer demand is aggregated onto the base bus (TBD with infrastructure PR)
    df = pd.concat(dfs, axis=1).T.groupby(level=0).sum().T

    return align_demand_to_snapshots(df, snapshots, format="%d. %b. %H:%M")


def set_load_sign(
    MM_data: pd.DataFrame,
    tables: list = ["power_generation", "hydrogen_supply"],
) -> pd.DataFrame:
    """
    Set negative sign for load values in market model data.

    The timeseries sheets report flags if a carrier is generation or load, with 
    load being reported with positive values.

    Parameters
    ----------
    MM_data : pd.DataFrame
        Market model data with columns 'flow', 'table', and 'value'.
    tables : list, default ["power_generation", "hydrogen_supply"]
        List of table names to apply the sign conversion to.

    Returns
    -------
    pd.DataFrame
        DataFrame with negated values for identified load carriers.
    """
    load_i = MM_data[
        (MM_data.flow == "Load") & MM_data.table.isin(tables)
    ].index
    MM_data.loc[load_i, "value"] *= -1
    return MM_data


def clean_h2_imports_for_benchmarking(
    crossborder_h2: pd.DataFrame,
    eu27: list[str],
) -> pd.DataFrame:
    """
    Extracts H2 imports from external nodes/buses for hydrogen supply benchmarking.

    Filters the hydrogen imports of the "Exchanges" sheet and maps their node to a benchmark carriers:
    - "Ammonia_": "ammonia imports"
    - Other import nodes (eg. DZ): "imports (renewable & low carbon)"

    The importing node is "bus1" and exporting node is"bus0" for shipped ammonia imports. 
    One exception is the IB_SKh2E-UA element where the reported value is negative and the 
    flow direction name convention reversed. Therefore, the sheet's zone is taken for bus0
    as this is the importing node.

    Parameters
    ----------
    eu27 : list[str]
        List of EU27 country codes (2-letter ISO).

    Returns
    -------
    pd.DataFrame
        DataFrame with columns [carrier, bus, unit, table, value] for each importing country and an EU27 aggregated row.
    """
    df = (
        crossborder_h2.query("carrier == 'H2_Imports'").assign(
            carrier=lambda x: np.where(
                x.element.str.contains("Ammonia"),
                "ammonia imports",
                "imports (renewable & low carbon)",
            ),
            unit="MWh",
            table="hydrogen_supply",
            value= lambda x: np.where(x.bus == x.bus1, x["sum"], -x["sum"])
        )
    )[["carrier", "bus", "value", "unit", "table"]]


    df_eu27 = (
        df[df.bus.map(extract_country).isin(eu27)]
        .groupby("carrier", as_index=False)
        .value.sum()
        .assign(bus="EU27", unit="MWh", table="hydrogen_supply")
    )
    return pd.concat([df, df_eu27], ignore_index=True)


def clean_crossborder_for_benchmarking(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean crossborder data for benchmarking purposes.
    """
    unit_map = df.loc["sum"].xs("unit", level="border")
    df = (
        df.loc["sum"]
        .T.reset_index()
        .query("border!='unit'")
        .rename(columns={"sum": "value"})
        .assign(
            unit=lambda x: x.carrier.map(unit_map),
            table=lambda x: np.where(
                x.carrier == "electricity",
                "crossborder_electricity",
                "crossborder_hydrogen",
            ),
            carrier=lambda x: np.where(x.carrier == "electricity", "AC", x.carrier),
        )
        .reset_index(drop=True)
    )

    df = (
        normalize_direction(df.set_index("border"), ["value"], buses_from_index=True)
        .reset_index()
        .drop(columns=["bus0", "bus1"])
    )
    return df


def assign_meta_data(df, planning_horizon, scenario):
    df["scenario"] = f"TYNDP {scenario}"
    df["year"] = planning_horizon
    df["source"] = "TYNDP 2026 Market Model Outputs"


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
    weather_scenario = get_weather_scenario(
        snakemake.params["weather_scenarios"], planning_horizon
    )

    # load carrier mapping
    mm_carrier_mapping = _load_mm_carrier_mapping(
        snakemake.input.carrier_mapping, options["tables"]
    )

    # EU27 country codes
    cc = coco.CountryConverter()
    eu27 = cc.EU27as("ISO2").ISO2.tolist()

    # TYNDP market model output files
    tyndp_output_file = sorted(
        Path(
            snakemake.input.tyndp_output_file,
            str(planning_horizon),
            f"Weather scenario {weather_scenario:03d}",
        ).glob("*_TimeSeriesDashboard_*.xlsx")
    )
    if not tyndp_output_file:
        raise FileNotFoundError(
            f"No TimeSeries Dashboard files for planning year {planning_horizon} "
            f"and weather scenario WS{weather_scenario:03d}."
        )

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
            mapping=mm_carrier_mapping,
        )

    MM_data = pd.concat(benchmarks).reset_index(drop=True)

    # set negative sign for loads
    MM_data = set_load_sign(MM_data)

    # load crossborder data
    logger.info(f"Processing tables of cross-border flows")
    crossborder = load_crossborder_sheet("Exchanges", tyndp_output_files)


    MM_data = pd.concat([MM_data, clean_crossborder_for_benchmarking(crossborder)])

    # concatenate crossborder H2 imports to MM data for benchmarking
    MM_data = pd.concat(
        [MM_data, clean_h2_imports_for_benchmarking(crossborder, eu27)],
        ignore_index=True,
    )

    # load h2 demand time series
    logger.info("Processing hourly H2 demand tables")
    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )
    h2_demand_ts = load_demand_ts(
        sheet_name="Hourly H2 Data",
        filepath=tyndp_output_file,
        snapshots=snapshots,
        carrier="H2",
    )
    elec_demand_ts = load_demand_ts(
        sheet_name="Hourly Market Data",
        filepath=tyndp_output_file,
        snapshots=snapshots,
        carrier="electricity",
    )

    # assign meta data
    assign_meta_data(MM_data, planning_horizon, scenario)
    assign_meta_data(crossborder, planning_horizon, scenario)
    assign_meta_data(h2_demand_ts, planning_horizon, scenario)

    # Save data
    MM_data.to_csv(snakemake.output.benchmarks, index=False)
    crossborder.to_csv(snakemake.output.crossborder)
    h2_demand_ts.to_csv(snakemake.output.h2_demand)
    elec_demand_ts.to_csv(snakemake.output.elec_demand)
