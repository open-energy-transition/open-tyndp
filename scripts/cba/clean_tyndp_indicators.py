# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Extract and clean 2024 TYNDP CBA indicators for transmission and storage projects.

Reads scenario sheets (tabs starting with 2030/2040), maps shortcut columns
using the Readme table (F6:I25), converts units to match Open-TYNDP
indicator units, and writes a long-format CSV for benchmarking.
"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

STAT_MAP = {
    "weighted avg": "mean",
    "avg": "mean",
    "max": "max",
    "min": "min",
    "low": "low",
    "mid": "central",
    "high": "high",
    "explicit": "explicit",
}

EURO_SYMBOLS = {
    "\u00e2\u0082\u00ac",
    "\u201a\u00c7\u00a8",
}

INDICATOR_MAP = {
    "B1": "B1_total_system_cost_change",
    "B3": "B3_res_generation_change_mwh",
    "B3a": "B3a_res_capacity_change_mw",
    "B4a": "B4a_nox",
    "B4b": "B4b_nh3",
    "B4c": "B4c_sox",
    "B4d": "B4d_pm25",
    "B4e": "B4e_pm10",
    "B4f": "B4f_nmvoc",
}

COMPOSITE_INDICATORS = {
    "B2a_co2_variation": {"keys": {"B2a", "B2b"}},
    "B2a_societal_cost_variation": {
        "keys": {
            "B2a_euro",
            "B2b_euro",
            "CO2_market_SEW",
            "CO2_network_SEW",
        }
    },
}

MODEL_UNITS = {
    "B1_total_system_cost_change": "EUR",
    "B2a_co2_variation": "t/year",
    "B2a_societal_cost_variation": "EUR/year",
    "B3_res_generation_change_mwh": "MWh/year",
    "B3a_res_capacity_change_mw": "MW",
    "B4a_nox": "kg/year",
    "B4b_nh3": "kg/year",
    "B4c_sox": "kg/year",
    "B4d_pm25": "kg/year",
    "B4e_pm10": "kg/year",
    "B4f_nmvoc": "kg/year",
}

UNIT_CONVERSION = {
    "Meuro/year": (1e6, "EUR/year"),
    "ktonnes/year": (1000.0, "t/year"),
    "GWh/year": (1000.0, "MWh/year"),
    "MWh/year": (1.0, "MWh/year"),
    "MW": (1.0, "MW"),
    "kg/year": (1.0, "kg/year"),
}


def normalize_unit(unit: str) -> str:
    unit_str = str(unit).strip()
    unit_str = unit_str.replace("\u20ac", "euro")
    unit_str = unit_str.replace(" ", "")
    return unit_str


def normalize_shortcut(shortcut: str) -> str:
    shortcut_str = str(shortcut).strip()
    shortcut_str = shortcut_str.replace("\u0394", "")
    shortcut_str = shortcut_str.replace("\u20ac", "euro")
    for bad in EURO_SYMBOLS:
        shortcut_str = shortcut_str.replace(bad, "euro")
    shortcut_str = shortcut_str.replace("monetized", "monetised")
    return shortcut_str


def normalize_indicator_text(text: str) -> str:
    value = str(text).strip()
    value = value.replace("\u0394", "")
    value = value.replace("\u20ac", "euro")
    for bad in EURO_SYMBOLS:
        value = value.replace(bad, "euro")
    value = value.replace("euroeuro", "euro")
    return value


def normalize_indicator_key(text: str) -> str:
    value = normalize_indicator_text(text)
    value = value.replace(" ", "")
    return value


def parse_project_id(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    extracted = series.astype(str).str.extract(r"(\d+)", expand=False)
    extracted = pd.to_numeric(extracted, errors="coerce")
    return numeric.fillna(extracted)


def convert_values(values: pd.Series, units) -> pd.Series:
    if isinstance(units, pd.Series):
        unit_norm = units.map(normalize_unit)
        factors = unit_norm.map(lambda u: UNIT_CONVERSION.get(u, (1.0, u))[0])
        return pd.to_numeric(values, errors="coerce") * factors
    unit_norm = normalize_unit(units)
    factor, _ = UNIT_CONVERSION.get(unit_norm, (1.0, unit_norm))
    return pd.to_numeric(values, errors="coerce") * factor


def load_readme_mapping(excel_path: Path) -> pd.DataFrame:
    mapping = pd.read_excel(
        excel_path, sheet_name="Readme", usecols="F:I", skiprows=5, nrows=25
    )
    mapping = mapping.dropna(subset=["Shortcut", "Indicator"]).copy()
    return mapping


def scenario_sheets(excel_path: Path) -> list[str]:
    xl = pd.ExcelFile(excel_path)
    return [
        name
        for name in xl.sheet_names
        if name[:4].isdigit() and name.startswith(("2030", "2040"))
    ]


def parse_scenario_name(sheet_name: str) -> tuple[str, int]:
    year = int(sheet_name[:4])
    return sheet_name, year


def extract_sheet(excel_path: Path, sheet_name: str) -> pd.DataFrame:
    df = pd.read_excel(excel_path, sheet_name=sheet_name, header=[0, 1])
    df = normalize_headers(df)

    project_cols = [c for c in df.columns if c[0] == "Project ID"]
    method_cols = [c for c in df.columns if c[0] == "Chosen Approach"]
    if not project_cols:
        logger.warning("Sheet %s has no Project ID column", sheet_name)
        return pd.DataFrame()

    project_col = project_cols[0]
    method_col = method_cols[0] if method_cols else None

    df[project_col] = parse_project_id(df[project_col])
    df = df.loc[~df[project_col].isna()].copy()
    if df.empty:
        return pd.DataFrame()

    df[project_col] = df[project_col].astype(int)

    indicator_cols = [
        c
        for c in df.columns
        if c[0] not in {"Project ID", "Project name", "Chosen Approach"}
    ]
    if not indicator_cols:
        logger.warning("No indicator columns found in %s", sheet_name)
        return pd.DataFrame()

    long = df[indicator_cols].copy()
    long.index = df[project_col]
    long = long.stack(level=[0, 1], future_stack=True).reset_index()
    long.columns = ["project_id", "shortcut", "stat", "value"]
    long["value"] = pd.to_numeric(long["value"], errors="coerce")

    stat = long["stat"].astype(str).str.strip().str.lower()
    long["subindex"] = stat.map(STAT_MAP).fillna(stat)
    long.loc[long["subindex"].isin(["nan", ""]), "subindex"] = None
    long = long.dropna(subset=["subindex"])

    if method_col is not None:
        method = df[[project_col, method_col]].copy()
        method.columns = ["project_id", "method"]
        method["method"] = (
            method["method"]
            .astype(str)
            .str.lower()
            .str.replace("light ", "", regex=False)
            .str.upper()
        )
        long = long.merge(method, on="project_id", how="left")
    else:
        long["method"] = None

    return long


def normalize_headers(df: pd.DataFrame) -> pd.DataFrame:
    level0 = []
    prev = None
    for col in df.columns:
        top = col[0]
        if pd.isna(top) or str(top).startswith("Unnamed"):
            top = prev
        level0.append(top)
        prev = top
    df.columns = pd.MultiIndex.from_tuples(
        list(zip(level0, df.columns.get_level_values(1)))
    )
    return df


def build_benchmark_rows(long: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    rows = []

    long["shortcut_norm"] = long["shortcut"].map(normalize_shortcut)
    mapping = mapping.copy()
    mapping["Shortcut_norm"] = mapping["Shortcut"].map(normalize_shortcut)
    merged = long.merge(
        mapping[["Shortcut_norm", "Indicator", "Unit", "Full name"]],
        left_on="shortcut_norm",
        right_on="Shortcut_norm",
        how="left",
    )
    merged = merged.rename(
        columns={
            "Indicator": "indicator_raw",
            "Unit": "unit_raw",
            "Full name": "indicator_name",
        }
    )
    merged["indicator_raw"] = merged["indicator_raw"].fillna(merged["shortcut_norm"])
    merged["indicator_raw"] = merged["indicator_raw"].map(normalize_indicator_text)
    merged["indicator_name"] = merged["indicator_name"].map(normalize_indicator_text)
    merged["indicator_key"] = merged["indicator_raw"].map(normalize_indicator_key)
    merged["indicator_mapped"] = merged["indicator_key"].map(INDICATOR_MAP)
    merged["indicator"] = merged["indicator_raw"]

    merged["unit_raw"] = merged["unit_raw"].map(normalize_unit)
    merged["value_raw"] = merged["value"]
    merged["value_converted"] = convert_values(merged["value"], merged["unit_raw"])
    merged["value"] = merged["value_raw"]
    merged["units"] = merged["unit_raw"]

    rows.append(
        merged[
            [
                "project_id",
                "method",
                "indicator",
                "indicator_raw",
                "indicator_mapped",
                "subindex",
                "units",
                "unit_raw",
                "value",
                "value_converted",
                "indicator_name",
            ]
        ]
    )

    for indicator, spec in COMPOSITE_INDICATORS.items():
        keys = list(spec["keys"])
        subset = merged[merged["indicator_key"].isin(keys)].copy()
        if subset.empty:
            continue
        pivot = (
            subset.pivot_table(
                index=["project_id", "method", "subindex"],
                columns="indicator_key",
                values="value",
                aggfunc="sum",
            )
            .fillna(0)
            .reset_index()
        )
        present = [key for key in keys if key in pivot.columns]
        if not present:
            continue
        pivot["value"] = pivot[present].sum(axis=1)
        pivot["indicator"] = indicator
        pivot["indicator_raw"] = None
        pivot["indicator_mapped"] = indicator
        pivot["indicator_name"] = None
        pivot["units"] = MODEL_UNITS.get(indicator, "")
        pivot["unit_raw"] = None
        pivot["value_converted"] = pivot["value"]
        rows.append(
            pivot[
                [
                    "project_id",
                    "method",
                    "indicator",
                    "indicator_raw",
                    "indicator_mapped",
                    "subindex",
                    "units",
                    "unit_raw",
                    "value",
                    "value_converted",
                    "indicator_name",
                ]
            ]
        )

    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_tyndp_indicators", run="NT", configfiles=["config/config.tyndp.yaml"]
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    def process_excel(excel_path: Path, project_type: str) -> pd.DataFrame:
        mapping = load_readme_mapping(excel_path)
        frames = []
        for sheet_name in scenario_sheets(excel_path):
            scenario, planning_horizon = parse_scenario_name(sheet_name)
            long = extract_sheet(excel_path, sheet_name)
            if long.empty:
                continue
            bench = build_benchmark_rows(long, mapping)
            if bench.empty:
                continue
            bench["scenario"] = scenario
            bench["planning_horizon"] = planning_horizon
            bench["type"] = project_type
            frames.append(bench)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    transmission_path = Path(snakemake.input.transmission)
    storage_path = Path(snakemake.input.storage)

    output_frames = [
        process_excel(transmission_path, "transmission"),
        process_excel(storage_path, "storage"),
    ]
    output = pd.concat([df for df in output_frames if not df.empty], ignore_index=True)

    if output.empty:
        output = pd.DataFrame(
            columns=[
                "project_id",
                "method",
                "indicator",
                "indicator_raw",
                "indicator_mapped",
                "indicator_name",
                "subindex",
                "units",
                "unit_raw",
                "value",
                "value_converted",
                "scenario",
                "planning_horizon",
                "type",
            ]
        )

    output.to_csv(snakemake.output.indicators, index=False)
