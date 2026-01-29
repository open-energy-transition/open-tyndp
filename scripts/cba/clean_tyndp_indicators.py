# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Extract and clean 2024 TYNDP CBA indicators for transmission and storage projects.

Reads scenario sheets (tabs starting with 2030/2040), maps shortcuts using the
table on the "Readme" tab, and outputs a long-format CSV.
"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

# Mapping of naming conventions of indicator results to standardized names
STAT_MAP = {
    "weighted avg": "mean",
    "avg": "mean",
    "mid": "central",
}

# Mapping of indicator shortcut names to standardized indicator names
# These should match the indicator names used in the Open-TYNDP CBA
INDICATOR_MAP = {
    "B1": "B1_total_system_cost_change",
    "B2a": "B2a_co2_variation",
    "B2a_euro": "B2a_societal_cost_variation",
    "B3": "B3_res_generation_change",
    "B3a": "B3a_res_capacity_change",
    "B4a": "B4a_nox",
    "B4b": "B4b_nh3",
    "B4c": "B4c_sox",
    "B4d": "B4d_pm25",
    "B4e": "B4e_pm10",
    "B4f": "B4f_nmvoc",
}


def normalize_text(
    value: str, *, drop_spaces: bool = False, monetised: bool = False
) -> str:
    """
    Normalize text for parsing the TYNDP Excel data.

    Strips whitespace, remove Delta symbol, replace Euro symbols with 'euro',
    standardizes spelling.
    """
    text = (
        str(value)
        .strip()
        .replace("\u0394", "")
        .replace("\u20ac", "euro")
        .replace("\u00e2\u0082\u00ac", "euro")
        .replace("\u201a\u00c7\u00a8", "euro")
    )
    if monetised:
        text = text.replace("monetized", "monetised")
    if drop_spaces:
        text = text.replace(" ", "")
    return text


def normalize_unit(unit: str) -> str:
    """Normalize unit strings."""
    return normalize_text(unit, drop_spaces=True)


def normalize_shortcut(shortcut: str) -> str:
    """Normalize Readme shortcut text."""
    text = normalize_text(shortcut, monetised=True)
    text = text.replace("CO2_market_SEW", "CO2_market_monetised")
    text = text.replace("CO2_network_SEW", "CO2_network_monetised")
    return text


def normalize_indicator_text(text: str) -> str:
    """Normalize indicator text."""
    return normalize_text(text)


def normalize_indicator_key(text: str) -> str:
    """Normalize indicator key."""
    return normalize_text(text, drop_spaces=True)


def parse_project_id(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    extracted = series.astype(str).str.extract(r"(\d+)", expand=False)
    extracted = pd.to_numeric(extracted, errors="coerce")
    return numeric.fillna(extracted)


def extract_readme_table(excel_path: Path) -> pd.DataFrame:
    """
    Load and normalize the Readme mapping table (F6:I25) from the TYNDP workbook.

    The Readme table is a table of indicator names, descriptions, shortcuts (shorthand name),
    and units. Returns a normalized DataFrame for export as CSV.
    """
    mapping = pd.read_excel(
        excel_path,
        sheet_name="Readme",
        usecols="F:I",
        skiprows=5,
        nrows=20,
    )
    mapping.columns = [str(col).strip() for col in mapping.columns]
    required = {"Shortcut", "Full name", "Indicator", "Unit"}
    missing = required - set(mapping.columns)
    if missing:
        logger.warning("Readme mapping missing columns: %s", sorted(missing))
        return pd.DataFrame(columns=["shortcut", "full_name", "indicator", "unit"])
    cleaned = mapping.rename(
        columns={
            "Indicator": "indicator",
            "Shortcut": "shortcut",
            "Unit": "unit",
            "Full name": "full_name",
        }
    )
    cleaned["indicator"] = cleaned["indicator"].map(normalize_indicator_text)
    cleaned["shortcut"] = cleaned["shortcut"].map(normalize_shortcut)
    cleaned["unit"] = cleaned["unit"].map(normalize_unit)
    cleaned["full_name"] = cleaned["full_name"].map(normalize_indicator_text)
    return cleaned[["shortcut", "full_name", "indicator", "unit"]]


def get_scenarios_sheets(excel_path: Path) -> list[str]:
    """
    Automatically detect and list scenario sheets in the Excel file.

    This assumes scenarios are named starting with 2030 or 2040.
    """
    xl = pd.ExcelFile(excel_path)
    return [
        name
        for name in xl.sheet_names
        if name[:4].isdigit() and name.startswith(("2030", "2040"))
    ]


def parse_scenario_name(sheet_name: str) -> tuple[str, int]:
    """Parse scenario name into label and planning horizon."""
    year = int(sheet_name[:4])
    return sheet_name, year


def extract_sheet(excel_path: Path, sheet_name: str) -> pd.DataFrame:
    """Extract long-format indicator rows from a single scenario sheet."""
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
    """Fill forward top-level headers in a multi-index header row."""
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


def join_indicator_data(long: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    """
    Join Readme metadata onto indicator data.

    Some data processing and normalization is applied to match the indicators
    between the two tables.
    """
    long["shortcut_norm"] = long["shortcut"].map(normalize_shortcut)
    mapping = mapping.copy()
    merged = long.merge(
        mapping[["shortcut", "indicator", "unit", "full_name"]],
        left_on="shortcut_norm",
        right_on="shortcut",
        how="left",
    )
    merged = merged.rename(
        columns={
            "indicator": "indicator_raw",
            "unit": "unit_raw",
            "full_name": "indicator_name",
        }
    )
    merged["indicator_raw"] = merged["indicator_raw"].fillna(merged["shortcut_norm"])
    merged["indicator_raw"] = merged["indicator_raw"].map(normalize_indicator_text)
    merged["indicator_name"] = merged["indicator_name"].map(normalize_indicator_text)
    merged["indicator_key"] = merged["indicator_raw"].map(normalize_indicator_key)
    merged["indicator_mapped"] = merged["indicator_key"].map(INDICATOR_MAP)
    merged["indicator"] = merged["indicator_raw"]

    merged["unit_raw"] = merged["unit_raw"].map(normalize_unit)
    merged["units"] = merged["unit_raw"]

    cols = [
        "project_id",
        "method",
        "indicator",
        "indicator_raw",
        "indicator_mapped",
        "subindex",
        "units",
        "unit_raw",
        "value",
        "indicator_name",
    ]
    base = merged[cols]

    return base


def process_benchmark_excel(excel_path: Path, project_type: str) -> pd.DataFrame:
    """Process all scenario sheets for a given workbook and project type."""
    mapping = extract_readme_table(excel_path)
    frames = []
    for sheet_name in get_scenarios_sheets(excel_path):
        scenario, planning_horizon = parse_scenario_name(sheet_name)
        long = extract_sheet(excel_path, sheet_name)
        if long.empty:
            continue
        bench = join_indicator_data(long, mapping)
        if bench.empty:
            continue
        bench["scenario"] = scenario
        bench["planning_horizon"] = planning_horizon
        bench["type"] = project_type
        frames.append(bench)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_tyndp_indicators", run="NT", configfiles=["config/config.tyndp.yaml"]
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    transmission_path = Path(snakemake.input.transmission)
    storage_path = Path(snakemake.input.storage)

    output_frames = [
        process_benchmark_excel(transmission_path, "transmission"),
        process_benchmark_excel(storage_path, "storage"),
    ]
    output = pd.concat([df for df in output_frames if not df.empty], ignore_index=True)

    output.to_csv(snakemake.output.indicators, index=False)

    readme_table = extract_readme_table(transmission_path)
    readme_table.to_csv(snakemake.output.readme, index=False)
    logger.info("Saved TYNDP indicator metadata to %s", snakemake.output.readme)
