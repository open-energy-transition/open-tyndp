# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script cleans the TYNDP market model output data for benchmarking.

Reads TYNDP market model outputs from the TimeSeries Dashboard xlsx files, one per country,
per planning horizon and weather scenario
- "Installed Capacity" sheet for installed capacities per market zone
- Nodal sheets for aggregates of hourly generation, load, prices, curtailment and unserved energy
- "Exchanges" sheet for cross-border electricity and hydrogen flows

Each file is read once and the files are processed in parallel.

The benchmark tables are defined in ``LOOKUP_TABLES`` as ``{table: options}`` with:
- ``sheet`` (required): Dashboard sheets to read the table from.
- ``stats`` (required, except for "power_capacity"): Yearly aggregate to use, one of
  "min", "max", "avg" or "sum".
- ``category`` (optional): Dashboard categories (without unit) to keep. All categories
  are kept if omitted.
- ``carrier`` (optional): Benchmark carrier assigned to all values. If omitted, carriers
  are mapped with the ``mapping_col`` of the table in ``tyndp_technology_map.csv``.

Note: Currently, only NT scenario processing is supported.
"""

import logging
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import country_converter as coco
import numpy as np
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.build_tyndp_network import extract_country
from scripts.sb._helpers import (
    CAPACITY_SHEET,
    ELECTRICITY_SHEETS,
    EXCHANGES_SHEET,
    H2_SHEETS,
    HEAT_SHEETS,
    NROWS_HEADER,
    STATS,
    convert_units,
    format_bz_names,
    get_dashboard_files,
    get_wscenario,
    normalize_direction,
    parse_dashboard_header,
    parse_installed_capacity,
)

logger = logging.getLogger(__name__)


LOOKUP_TABLES: dict[str, dict] = {
    "power_capacity": {"sheet": [CAPACITY_SHEET]},
    "power_generation": {"sheet": ELECTRICITY_SHEETS, "stats": "sum"},
    # E-Market and Prosumer demand share one benchmark carrier in current build_statistics
    "electricity_demand": {
        "sheet": ELECTRICITY_SHEETS,
        "category": [
            "Native Demand",
            "Fixed Demand",
            "CH4 Heat Pump (load)",
            "H2 Heat Pump (load)",
        ],
        "stats": "sum",
    },
    "electricity_prosumer_demand": {
        "sheet": ["Prosumer"],
        "category": ["Native Demand", "Fixed Demand"],
        "stats": "sum",
    },
    "hydrogen_demand": {
        "sheet": H2_SHEETS,
        "category": [
            "Native Demand",
            "H2 to H2 power plants",
            "H2 to eLiquids",
            "H2 to sng",
            "H2 Boiler (load)",
        ],
        "stats": "sum",
    },
    "hydrogen_supply": {"sheet": H2_SHEETS, "stats": "sum"},
    "hybrid_heating": {
        "sheet": HEAT_SHEETS,
        "category": ["CH4 Heat Pump", "CH4 Boiler", "H2 Heat Pump", "H2 Boiler"],
        "stats": "sum",
    },
    "ev_charging": {
        "sheet": ["E-Market", "Prosumer"],
        "category": ["EV charge (load)", "EV discharge (gen.)"],
        "stats": "sum",
    },
    # TODO: For shedding hours, TYNDP 2026 does not give a direct value, so it can be derived
    # where "Energy Not Served" is above 0.
    # "electricity_demand_shedding_hours": [
    # "Yearly Outputs",
    # "Loss of load expectation [hour]  ",
    # ],  # includes white space
    # "hydrogen_demand_shedding_hours": [
    # "Yearly H2 Outputs",
    # "Loss of H2 load expectation [hour]  ",
    # ],  # includes white space
    # prices
    "electricity_price": {
        "sheet": ["E-Market"],
        "category": ["Marginal Cost"],
        "stats": "avg",
        "carrier": "AC",
    },
    # "electricity_price_excl_shed": [
    # "Yearly Outputs",
    # "Marginal Cost Yearly Average (excl. 3 000 €/MWh) [€]",
    # ],
    "hydrogen_price": {
        "sheet": ["H2 Zone 2"],
        "category": ["Marginal Cost"],
        "stats": "avg",
        "carrier": "H2",
    },
    # "hydrogen_price_excl_shed": [
    # "Yearly H2 Outputs",
    # "Marginal Cost Yearly Average (excl. 3 000 €/MWhH2) [€/MWhH2]",
    # ],
}

# look up dictionary for crossborder exchanges
CROSS_BORDER_DICT: dict[str, str] = {
    "E-Market Exchanges": "electricity",
    "H2 Zone 2 Exchanges": "H2",
    "H2 Zone 2 to H2 Zone 1": "H2",
    "H2 Imports to H2 Zone 2": "H2_imports",
}


def _load_dashboard_carrier_mapping(
    carrier_mapping_fn: str, tables: dict
) -> dict[str, dict]:
    """
    Load mapping from TYNDP dashboard carrier names to benchmark carrier names.
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
            tech_map[["tyndp_dashboard_carrier", col]]
            .dropna(subset=["tyndp_dashboard_carrier", col])
            .set_index("tyndp_dashboard_carrier")[col]
            .to_dict()
        )

    return output_map


def rename_prosumer_nodes(buses: pd.Series) -> pd.Series:
    """
    Rename prosumer and hybrid heating nodes to match Open-TYNDP bus names
    """
    return format_bz_names(buses).str.replace(r"(RETE| HCH4| HH2)$", "", regex=True)


def parse_dashboard_sheet(df: pd.DataFrame, sheet_name: str) -> pd.DataFrame:
    """
    Parse the yearly aggregates of a sheet from a TYNDP TimeSeries Dashboard output file.

    See ``parse_dashboard_header`` for the structure of the sheet.

    Parameters
    ----------
    df : pd.DataFrame
        Raw sheet read without header.
    sheet_name : str
        Name of the Excel sheet.

    Returns
    -------
    pd.DataFrame
        Long format with columns [bus, element, flow, carrier, unit, stats,
        value, sheet]
    """
    df = (
        parse_dashboard_header(df)
        .melt(
            id_vars=["bus", "element", "flow", "carrier", "unit"],
            value_vars=STATS,
            var_name="stats",
        )
        .assign(sheet=sheet_name)
    )
    sum_unit = "TJ" if sheet_name in HEAT_SHEETS else "GWh"
    df["unit"] = df.unit.where(df.stats != "sum", sum_unit)
    df["value"] = pd.to_numeric(df["value"], errors="coerce")

    return df


def read_dashboard(
    filepath: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read all required sheets of a TYNDP TimeSeries Dashboard output file at once.

    Parameters
    ----------
    filepath : str or Path
        Path to the Excel file

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        Yearly aggregates of all timeseries sheets and installed capacities.
    """
    sheet_names = ELECTRICITY_SHEETS + H2_SHEETS + HEAT_SHEETS + [EXCHANGES_SHEET]

    # Some nodes are missing sheets
    with pd.ExcelFile(filepath, engine="calamine") as file:
        sheets = pd.read_excel(
            file,
            sheet_name=[s for s in sheet_names if s in file.sheet_names],
            header=None,
            nrows=NROWS_HEADER,
        )
        capacity = pd.read_excel(file, sheet_name=CAPACITY_SHEET, header=None)

    stats = pd.concat([parse_dashboard_sheet(df, s) for s, df in sheets.items()])

    return stats, parse_installed_capacity(capacity)


def drop_single_zone_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Drop the cross-border flows reported by a country as a single zone.

    In the TYNDP 2026 dashboard outputs for 2050, Italy reports its exchanges
    as a single zone `IT`, e.g. `AT00-IT`, while its neighbours report the
    same flows with the Italian zones, e.g. `AT00-ITN1`, which match the
    model borders.

    Parameters
    ----------
    df : pd.DataFrame
        Cross-border flows with the reporting zone as `bus` and the border as
        `border`.

    Returns
    -------
    pd.DataFrame
        Cross-border flows without the single zone duplicates.
    """
    single_zones = {
        bus
        for bus in df.loc[df.carrier.eq("electricity"), "bus"].unique()
        if bus == extract_country(bus)
    }
    ends = df.border.str.split("-")
    return df[~ends.map(lambda x: bool(single_zones.intersection(x)))]


def load_crossborder(
    df: pd.DataFrame,
    stats: list[str] = STATS,
) -> pd.DataFrame:
    """
    Load the cross-border flows from the "Exchanges" sheets of the TYNDP TimeSeries Dashboard.

    In TYNDP 2026, both electricity and hydrogen carriers share the "Exchanges" sheet.

    Parameters
    ----------
    df : pd.DataFrame
        Yearly aggregates of all timeseries sheets in long format.
    stats : list[str], optional
        Yearly aggregates to extract

    Returns
    -------
    pd.DataFrame
        DataFrame with normalized cross-border flow data.
    """
    df = (
        df.query("sheet == @EXCHANGES_SHEET and stats in @stats")
        .assign(
            carrier=lambda x: x.carrier.map(CROSS_BORDER_DICT),
            border=lambda x: format_bz_names(x.element),
        )
        .dropna(subset=["carrier"])
        .pipe(drop_single_zone_duplicates)
        .drop_duplicates(["border", "stats"])
        .pivot(index=["carrier", "bus", "border"], columns="stats", values="value")
        .reset_index(["carrier", "bus"])
    )

    # normalize direction
    df[["bus0", "bus1"]] = df.index.to_series().str.split("-", expand=True).values
    df = normalize_direction(df, cols=stats).set_index("border")
    mask = df["min"] > df["max"]
    df.loc[mask, ["min", "max"]] = df.loc[mask, ["max", "min"]].values

    # convert units
    df["sum"] = df["sum"].mul(1e3)

    return df.sort_index()


def load_dashboard_table(
    table_name: str,
    stats: pd.DataFrame,
    capacity: pd.DataFrame,
    countries: list[str],
    eu27: list,
    mapping: dict[str, dict[str, str]],
) -> pd.DataFrame:
    """
    Build a benchmarking table from TYNDP 2026 TimeSeries Dashboard data.

    Parameters
    ----------
    table_name : str
        Name of the table from LOOKUP_TABLES (e.g., "power_capacity").
    stats : pd.DataFrame
        Yearly aggregates of all timeseries sheets in long format.
    capacity : pd.DataFrame
        Installed capacities in long format.
    countries : list[str]
        List of modelled countries.
    eu27 : list
        List of EU27 country codes.
    mapping : dict[str, dict[str, str]]
        Carrier mapping from dashboard carrier names to benchmarking carrier names per table.

    Returns
    -------
    pd.DataFrame
        Dashboard data in long format (incl. EU27).
    """
    opt = LOOKUP_TABLES[table_name]

    if table_name == "power_capacity":
        df = capacity.copy()
    else:
        df = stats[stats.sheet.isin(opt["sheet"]) & (stats.stats == opt["stats"])]
        if categories := opt.get("category"):
            df = df[df.carrier.isin(categories)]

    # Only include mapped carriers
    if carrier := opt.get("carrier"):
        df = df.assign(carrier=carrier)
    else:
        carriers = df.carrier.map(mapping[table_name])
        if unmapped := sorted(df.carrier[carriers.isna()].unique()):
            logger.warning(
                f"No carrier mappings for table '{table_name}' for: {unmapped}. They will be excluded from the benchmarking for this table."
            )
        df = df.assign(carrier=carriers).dropna(subset=["carrier"])

    df = set_load_sign(df, table_name)

    # Rename and filter column names (buses)
    df = df.assign(bus=lambda x: rename_prosumer_nodes(x.bus.astype(str)))
    op = "sum" if "price" not in table_name else "mean"
    df_nodal = df.groupby(["bus", "carrier", "unit"], as_index=False).value.agg(op)
    df_nodal = df_nodal[df_nodal.bus.map(extract_country).isin(countries)]

    # Add EU27 / Pan-EU load-weighted average for prices
    if "price" in table_name:
        df_eu = df_nodal
        bus_name = "Pan-EU"
        weights = (
            load_dashboard_table(
                table_name=f"{table_name.split('_')[0]}_demand",
                stats=stats,
                capacity=capacity,
                countries=countries,
                eu27=eu27,
                mapping=mapping,
            )
            .groupby("bus")
            .value.sum()
            .reindex(df_eu.bus.unique(), fill_value=0)
        )
        normalizer = weights.sum()
    else:
        df_eu = df_nodal[df_nodal.bus.map(extract_country).isin(eu27)]
        bus_name = "EU27"
        weights = pd.Series(1.0, index=df_eu.bus.unique())
        normalizer = 1

    df_eu = (
        df_eu.assign(value=lambda x: x.bus.map(weights) * x.value)
        .groupby(by=["carrier", "unit"])
        .value.sum()
        .div(normalizer)
        .reset_index()
        .assign(bus=bus_name)
    )
    df = pd.concat([df_nodal, df_eu])

    df["table"] = table_name
    if "price" not in table_name:
        df = convert_units(df)

    return df


def set_load_sign(
    df: pd.DataFrame,
    table_name: str,
    tables: list = ["power_generation", "hydrogen_supply"],
) -> pd.DataFrame:
    """
    Set negative sign for load values in market model data.

    The timeseries sheets report flags if a carrier is generation or load, with
    load being reported with positive values.

    Parameters
    ----------
    df : pd.DataFrame
        Market model data with columns 'flow' and 'value'.
    table_name : str
        Name of the table from LOOKUP_TABLES.
    tables : list, default ["power_generation", "hydrogen_supply"]
        List of table names to apply the sign conversion to.

    Returns
    -------
    pd.DataFrame
        DataFrame with negated values for identified load carriers.
    """
    if table_name not in tables:
        return df

    return df.assign(value=df.value.where(df.flow != "Load", -df.value))


def clean_crossborder_for_benchmarking(
    df: pd.DataFrame, eu27: list[str]
) -> pd.DataFrame:
    """
    Clean crossborder data for benchmarking purposes.
    """
    tables = {
        "electricity": "crossborder_electricity",
        "H2": "crossborder_hydrogen",
        "H2_imports": "crossborder_hydrogen",
    }

    df = df.reset_index().rename(columns={"sum": "value"})

    flows = df[df.carrier.isin(tables)].assign(
        table=lambda x: x.carrier.map(tables),
        carrier=lambda x: x.carrier.replace({"electricity": "AC", "H2_imports": "H2"}),
    )[["border", "carrier", "value", "table"]]

    imports = df[df.carrier == "H2_imports"].assign(
        carrier=lambda x: np.where(
            x.border.str.contains("Ammonia"),
            "ammonia imports",
            "imports (renewable & low carbon)",
        ),
        value=lambda x: np.where(x.bus == x.bus1, x.value, -x.value),
        table="hydrogen_supply",
    )[["bus", "carrier", "value", "table"]]

    imports_eu27 = (
        imports[imports.bus.map(extract_country).isin(eu27)]
        .groupby(["carrier", "table"], as_index=False)
        .value.sum()
        .assign(bus="EU27")
    )

    df = pd.concat([flows, imports, imports_eu27], ignore_index=True)
    df["unit"] = "MWh"

    return df


def assign_meta_data(df, planning_horizon, scenario):
    df["scenario"] = f"TYNDP {scenario}"
    df["year"] = planning_horizon
    df["source"] = "TYNDP 2026 Dashboard Outputs"


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_tyndp_dashboard_benchmark",
            horizon="2040",
            scenario="NT",
            configfiles="config/test/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    options = snakemake.params["benchmarking"]
    scenario = snakemake.params["scenario"]
    planning_horizon = int(snakemake.wildcards.horizon)
    countries = snakemake.params["countries"]

    # currently only implemented for NT
    if scenario != "NT":
        logger.warning(
            "Processing of TYNDP dashboard files currently only implemented for NT scenario. Exporting empty Data Frames"
        )
        for fn in snakemake.output:
            pd.DataFrame().to_csv(fn)
        sys.exit(0)

    wscenario = get_wscenario(snakemake.params["wscenarios"], planning_horizon)

    # Logs which weather scenario and planning year is being processed
    logger.info(
        f"Processing planning year {planning_horizon}, weather scenario WS{wscenario:03d}"
    )

    # load carrier mapping
    carrier_mapping = _load_dashboard_carrier_mapping(
        snakemake.input.carrier_mapping, options["tables"]
    )

    # EU27 country codes
    cc = coco.CountryConverter()
    eu27 = cc.EU27as("ISO2").ISO2.tolist()

    # TYNDP TimeSeries Dashboard files
    dashboard_files = get_dashboard_files(
        snakemake.input.dashboard_dir, planning_horizon, wscenario
    )

    logger.info(f"Reading {len(dashboard_files)} TimeSeries Dashboard files")
    with ProcessPoolExecutor(max_workers=snakemake.threads) as executor:
        stats, capacity = zip(*executor.map(read_dashboard, dashboard_files))
    stats = pd.concat(stats, ignore_index=True)
    capacity = pd.concat(capacity, ignore_index=True)

    # Tables for which TYNDP dashboard files provide data
    tables_to_process = [
        t for t in LOOKUP_TABLES.keys() if t in options["tables"].keys()
    ]

    logger.info(f"Processing tables: {', '.join(tables_to_process)}")

    dashboard_data = pd.concat(
        [
            load_dashboard_table(
                table_name=table,
                stats=stats,
                capacity=capacity,
                countries=countries,
                eu27=eu27,
                mapping=carrier_mapping,
            )
            for table in tables_to_process
        ],
        ignore_index=True,
    )

    # load crossborder data
    logger.info("Processing tables of cross-border flows")
    crossborder = load_crossborder(stats)

    # concatenate crossborder flows and imports to dashboard_data for benchmarking
    dashboard_data = pd.concat(
        [dashboard_data, clean_crossborder_for_benchmarking(crossborder, eu27)],
        ignore_index=True,
    )

    # assign meta data
    assign_meta_data(dashboard_data, planning_horizon, scenario)

    # Save data
    dashboard_data.to_csv(snakemake.output.benchmarks, index=False, float_format="%.2f")
