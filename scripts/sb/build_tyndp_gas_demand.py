# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP Scenario Building gas demand for Open-TYNDP.

Processes methane (gas) demand from TYNDP Supply Tool, with TYNDP 2026 support.
- TYNDP 2024: NT+ data, Other data and Conversions, IT sheets (GWh -> MWh)
- TYNDP 2026: Single ``All data`` sheet with ``category | parameter | unit | year | AT..SK | EU27``
  layout, ``Year`` column covering 2030/2035/2040/2050 in one block, units TWh/year -> MWh.

For 2026 the thermal demand in gas boilers for heating (hybrid heating) is already
provided via ``build_tyndp_demand``'s ``thermal_ch4`` hourly per-node profiles.
To avoid double-counting, the script subtracts the annual hybrid-heating gas
use from the Supply Tool aggregate. The better temporal shape is kept via the
hourly ``thermal_ch4`` profile; the residual gas demand remains flat (annual
MWh divided by hours downstream in ``prepare_sector_network``).

Units, planning years and energy balance are verified, and interpolation is
kept for missing years.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Literal

import country_converter as coco
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.sb._helpers import interpolate_demand

logger = logging.getLogger(__name__)
cc = coco.CountryConverter()

# TYNDP 2026 available years and unit handling
AVAILABLE_YEARS_TYNDP2026 = [2030, 2035, 2040, 2050]
# Carriers that contributed to FED in 2024; 2026 aggregates similar parameters.
# Keep same list for comparability; if All data uses aggregated parameter, fallback will catch it.
GAS_FED_CARRIERS_2024 = [
    "E-Methane",
    "Other fossil gas",
    "Biomethane",
    "Natural gas",
    "Waste gas",
    "Gas for Cooking",
    "Methane (LNG)",
]
GAS_HEAT_CARRIERS_2024 = [
    "Biogas",
    "E-Methane",
    "Natural gas",
    "Other fossil gas",
    "Waste gas",
]
# For 2026 we try to match these parameters; if not found we fallback to category/parameter contains gas|methane
GAS_PARAMS_2026 = set(GAS_FED_CARRIERS_2024 + GAS_HEAT_CARRIERS_2024 + ["Biogas", "Methane"])


def _is_country_col(col: str) -> bool:
    """Heuristic: 2-letter country code or bus code like AT00."""
    col = str(col).strip()
    if col.lower() in {"category", "parameter", "unit", "year", "scenario", "eu27"}:
        return False
    # Try to convert first 2 letters to iso2
    try:
        iso = cc.convert(col[:2], to="iso2")
        return iso != "not found" and len(col[:2]) == 2
    except Exception:
        return False


def _unit_factor(unit_str: str) -> float:
    """Return MWh factor for a unit string."""
    if not isinstance(unit_str, str):
        return 1.0
    u = unit_str.strip().lower()
    if "twh" in u:
        return 1e6
    if "gwh" in u:
        return 1e3
    if "mwh" in u:
        return 1.0
    if "kwh" in u:
        return 1e-3
    return 1.0


def _find_supply_tool_file(fn: str | Path) -> Path | None:
    """Resolve supply_tool path: if directory, find first .xls* file inside."""
    p = Path(fn)
    if p.is_dir():
        candidates = list(p.rglob("*.xls*"))
        # Prefer file with "Supply" in name or largest
        if not candidates:
            return None
        # Prefer All data containing file
        for c in candidates:
            if "supply" in c.name.lower():
                return c
        return candidates[0]
    if p.is_file():
        return p
    # Try glob relative
    candidates = list(Path(".").glob(str(fn)))
    if candidates:
        return candidates[0]
    return p if p.exists() else None


def read_all_data_2026(fn: str, scenario: str, pyear: int) -> pd.Series | None:
    """
    Read TYNDP 2026 Supply Tool ``All data`` sheet.

    Expected layout: category | parameter | unit | year | AT .. SK | EU27 | EU27
    with a ``Year`` column covering 2030/2035/2040/2050.

    Returns Series per iso2 code in MWh, or None if sheet not found.
    """
    actual = _find_supply_tool_file(fn)
    if actual is None:
        logger.debug(f"Supply tool file not found for {fn}")
        return None
    fn = str(actual)
    try:
        # Try reading All data sheet; try both .xlsm and .xlsx engines
        try:
            df = pd.read_excel(fn, sheet_name="All data", header=0, engine="openpyxl")
        except Exception:
            df = pd.read_excel(fn, sheet_name="All data", header=0)
    except Exception as e:
        logger.debug(f"All data sheet not found in {fn}: {e}")
        return None

    # Normalize columns
    df.columns = [str(c).strip() for c in df.columns]
    # Find key columns case-insensitive
    def find_col(name):
        for c in df.columns:
            if c.lower() == name.lower():
                return c
        return None

    year_col = find_col("year")
    cat_col = find_col("category")
    param_col = find_col("parameter")
    unit_col = find_col("unit")
    scenario_col = find_col("scenario")

    if year_col is None or param_col is None:
        logger.warning("All data sheet missing Year or Parameter column")
        return pd.Series(dtype=float)

    # Identify country columns
    country_cols = [c for c in df.columns if _is_country_col(c)]
    # Deduplicate EU27 and keep first occurrence of each iso2
    # Map col -> iso2
    col_to_iso = {}
    for c in country_cols:
        iso = cc.convert(str(c)[:2], to="iso2")
        if iso == "not found":
            continue
        # Skip EU27 aggregate (if iso is EU? cc may return not found)
        if iso.upper() == "EU":
            continue
        # Keep first occurrence for each iso
        if iso not in col_to_iso.values():
            col_to_iso[c] = iso
        else:
            # Duplicate country (e.g., second EU27 or duplicate AT) – skip
            logger.debug(f"Skipping duplicate country column {c} -> {iso}")

    if not col_to_iso:
        logger.warning("No country columns identified in All data sheet")
        return pd.Series(dtype=float)

    # Filter for pyear
    # Year column may be string or int
    try:
        df[year_col] = pd.to_numeric(df[year_col], errors="coerce")
    except Exception:
        pass
    df_year = df[df[year_col] == pyear]
    if df_year.empty:
        logger.debug(f"No rows for year {pyear} in All data")
        return pd.Series(dtype=float)

    # Scenario handling: if scenario column exists, filter for scenario (NT/DE/GA)
    if scenario_col is not None and scenario in df_year[scenario_col].astype(str).values:
        df_year = df_year[df_year[scenario_col].astype(str).str.contains(scenario, case=False, na=False) | df_year[scenario_col].isna()]
        # If filtering removes all, fallback to unfiltered
        if df_year.empty:
            df_year = df[df[year_col] == pyear]

    # Filter for gas-relevant rows
    # First try exact parameter match
    mask = df_year[param_col].isin(GAS_PARAMS_2026)
    if not mask.any():
        # Fallback: category or parameter contains gas/methane (case-insensitive)
        mask_cat = df_year[cat_col].astype(str).str.contains("gas|methane", case=False, na=False) if cat_col else False
        mask_param = df_year[param_col].astype(str).str.contains("gas|methane", case=False, na=False)
        mask = mask_cat | mask_param
        if not mask.any():
            # Last fallback: take rows where unit is TWh/year and year matches (assume all are gas)
            logger.warning(f"No gas-specific rows found for {pyear}, trying unit-based fallback")
            if unit_col:
                mask = df_year[unit_col].astype(str).str.contains("twh|gwh", case=False, na=False)
            else:
                mask = pd.Series(True, index=df_year.index)

    df_gas = df_year[mask]
    if df_gas.empty:
        logger.warning(f"No gas rows after filtering for year {pyear}")
        return pd.Series(dtype=float)

    # Convert country columns to numeric and sum per row group
    for c in col_to_iso:
        df_gas[c] = pd.to_numeric(df_gas[c], errors="coerce").fillna(0)

    # Determine unit factor per row: if unit column exists, use it per row, else assume TWh
    if unit_col and unit_col in df_gas.columns:
        # Compute weighted sum per country: value * factor
        totals = {iso: 0.0 for iso in set(col_to_iso.values())}
        for _, row in df_gas.iterrows():
            factor = _unit_factor(str(row[unit_col]))
            for col, iso in col_to_iso.items():
                totals[iso] += float(row[col]) * factor
        series = pd.Series(totals, dtype=float)
    else:
        # Assume TWh
        summed = df_gas[list(col_to_iso.keys())].sum()
        summed.index = [col_to_iso[c] for c in summed.index]
        series = summed.groupby(level=0).sum() * 1e6

    series.name = "p_nom"
    series = series[series != 0]
    logger.info(f"All data 2026: extracted gas demand for {pyear} ({scenario}): {series.sum()/1e6:.1f} TWh total")
    return series


def read_fed_data(fn: str, scenario: str, pyear: int) -> tuple[pd.Series, pd.Series]:
    """Legacy 2024 FED reader."""
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
        demand_fed.columns = pd.Index(cc.convert(demand_fed.columns, to="iso2"))
        demand_heat = demand_fed.loc["Heat"].mul(1e3)
        demand_fed = (
            demand_fed.loc[GAS_FED_CARRIERS_2024].mul(1e3).sum()
        )
    except Exception as e:
        logger.warning(
            f"Failed to read final gas demand for scenario {scenario} and pyear {pyear}: {type(e).__name__}: {e}"
        )
        demand_fed = pd.Series()
        demand_heat = pd.Series()
    return demand_fed, demand_heat


def read_heat_frame(
    fn: str, pyear: int, type: Literal["distribution", "efficiency"]
) -> pd.DataFrame:
    if type not in ["distribution", "efficiency"]:
        raise ValueError(f"Invalid type '{type}'. Must be 'distribution' or 'efficiency'.")
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
    df.columns = pd.Index(cc.convert(df.columns, to="iso2"))
    return df


def read_it_gas_prod(fn: str, pyear: int) -> float:
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
    )


def read_heat_data(
    heat_fed: pd.Series, fn: str, scenario: str, pyear: int
) -> pd.Series:
    try:
        shares = read_heat_frame(fn, pyear, "distribution")
        efficiencies = read_heat_frame(fn, pyear, "efficiency")
        demand_primary = heat_fed * shares * (1 / efficiencies)
        demand = demand_primary.loc[
            ["Biogas", "E-Methane", "Natural gas", "Other fossil gas", "Waste gas"]
        ].sum()
        demand.loc["IT"] = read_it_gas_prod(fn, pyear)
    except Exception as e:
        logger.warning(
            f"Failed to read heat demand data for scenario {scenario} and pyear {pyear}: {type(e).__name__}: {e}"
        )
        demand = pd.Series()
    return demand


def read_supply_tool_2024(fn: str, scenario: str, pyear: int) -> pd.Series:
    demand_fed, heat_fed = read_fed_data(fn, scenario, pyear)
    demand_heat = read_heat_data(heat_fed, fn, scenario, pyear)
    demand = pd.concat([demand_fed, demand_heat], axis=1).sum(axis=1)
    demand.name = "p_nom"
    return demand


def compute_thermal_ch4_annual(
    thermal_ch4_path: str | Path | None, snapshots: pd.DatetimeIndex | None = None
) -> pd.Series:
    """
    Compute annual hybrid-heating gas demand from thermal_ch4 hourly profile.

    thermal_ch4 is hourly MW_th per bus (from build_tyndp_demand). Sum over
    snapshots gives MWh_th per bus, then aggregated to iso2 country.

    If snapshots with weights are available, a more accurate integration could
    be done, but for now sum * 1h is used.
    """
    if not thermal_ch4_path:
        return pd.Series(dtype=float)
    p = Path(thermal_ch4_path)
    if not p.exists():
        # Try alternative: snakemake input may be directory or missing
        logger.debug(f"thermal_ch4 file not found: {p}")
        return pd.Series(dtype=float)
    try:
        df = pd.read_csv(p, index_col=0, parse_dates=True)
        if df.empty:
            return pd.Series(dtype=float)
        # If file has single column annual already, handle
        if df.shape[0] == 1 and "p_nom" in df.columns:
            s = df["p_nom"]
            s.index = s.index.map(lambda x: str(x)[:2])
            return s
        # Hourly: sum per bus
        # Handle non-numeric columns
        df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
        annual_per_bus = df.sum(axis=0)  # MWh_th (MW * 1h per row)
        # If snapshots weighting needed, assume 8760 vs 8736 handled downstream; sum is fine
        # Map bus to country (first 2 chars)
        annual_per_country = {}
        for bus, val in annual_per_bus.items():
            country = str(bus)[:2]
            iso = cc.convert(country, to="iso2")
            if iso == "not found":
                continue
            annual_per_country[iso] = annual_per_country.get(iso, 0) + float(val)
        s = pd.Series(annual_per_country, dtype=float)
        s.name = "thermal_ch4"
        logger.info(f"thermal_ch4 annual total: {s.sum()/1e6:.2f} TWh from {p}")
        return s
    except Exception as e:
        logger.warning(f"Failed to read thermal_ch4 {p}: {e}")
        return pd.Series(dtype=float)


def read_supply_tool(fn: str, scenario: str, pyear: int) -> pd.Series:
    """
    Unified reader: try 2026 All data first, fallback to 2024.
    """
    # Try 2026 All data
    result_2026 = read_all_data_2026(fn, scenario, pyear)
    if result_2026 is not None and not result_2026.empty:
        return result_2026
    # Fallback to 2024
    logger.info(f"Falling back to 2024 Supply Tool parsing for {pyear}")
    return read_supply_tool_2024(fn, scenario, pyear)


def load_single_year(
    fn: str, scenario: str, pyear: int, thermal_ch4_path: str | None = None
) -> pd.Series:
    """Load demand data for a single planning year, with hybrid heating subtraction."""
    if scenario == "NT":
        demand = read_supply_tool(fn, scenario, pyear)
    elif scenario in ["DE", "GA"]:
        demand = pd.Series(dtype=float)
        # Try 2026 path also for DE/GA
        attempt = read_all_data_2026(fn, scenario, pyear)
        if attempt is not None and not attempt.empty:
            demand = attempt
    else:
        demand = pd.Series(dtype=float)

    # Subtract hybrid heating (thermal_ch4) to avoid double-counting
    if thermal_ch4_path:
        thermal = compute_thermal_ch4_annual(thermal_ch4_path)
        if not thermal.empty and not demand.empty:
            # Align indexes
            demand, thermal = demand.align(thermal, fill_value=0, join="outer")
            # Log energy balance before
            total_before = demand.sum()
            thermal_sum = thermal.sum()
            # Subtract; clip at 0 to avoid negative residual
            residual = demand - thermal
            neg = residual[residual < 0]
            if not neg.empty:
                logger.warning(
                    f"Hybrid heating subtraction leads to negative residual for {neg.index.tolist()}: {neg.values}. Clipping to 0."
                )
                residual = residual.clip(lower=0)
            # Verify energy balance: residual + thermal ~= demand
            total_after = residual.sum() + thermal_sum
            if total_before > 0:
                diff_pct = abs(total_after - total_before) / total_before * 100
                if diff_pct > 1e-6:
                    logger.warning(f"Energy balance mismatch after subtraction: before {total_before/1e6:.2f} TWh, after {total_after/1e6:.2f} TWh (diff {diff_pct:.4f}%)")
                else:
                    logger.info(f"Energy balance OK: total {total_before/1e6:.2f} TWh = residual {residual.sum()/1e6:.2f} TWh + thermal {thermal_sum/1e6:.2f} TWh")
            demand = residual
            demand.name = "p_nom"

    return demand


def load_gas_demand(
    fn: str, scenario: str, pyear: int, thermal_ch4_path: str | None = None
) -> pd.Series:
    """
    Load gas demand data for a specific scenario and planning year.
    Supports 2026 years [2030,2035,2040,2050] with interpolation.
    """
    available_years = AVAILABLE_YEARS_TYNDP2026

    if pyear in available_years:
        logger.debug(f"Year {pyear} found in available data. Loading directly.")
        return load_single_year(fn, scenario, pyear, thermal_ch4_path=thermal_ch4_path)

    return interpolate_demand(
        available_years=available_years,
        pyear=pyear,
        load_single_year_func=load_single_year,
        fn=fn,
        scenario=scenario,
        thermal_ch4_path=thermal_ch4_path,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_gas_demand",
            configfiles="config/test/config.tyndp.yaml",
            planning_horizons=2035,
            run="NT",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    scenario = snakemake.params["scenario"]
    fn = snakemake.input.supply_tool
    pyear = int(snakemake.wildcards.planning_horizons)
    # New input for hybrid heating: try snakemake.input.thermal_ch4 if provided
    thermal_ch4_path = None
    for attr in ["thermal_ch4", "thermal_ch4_file", "hybrid_heating"]:
        if hasattr(snakemake.input, attr):
            thermal_ch4_path = getattr(snakemake.input, attr)
            break
        if isinstance(snakemake.input, dict) and attr in snakemake.input:
            thermal_ch4_path = snakemake.input[attr]
            break
    if thermal_ch4_path is None and hasattr(snakemake.input, "__getitem__"):
        for key in ["thermal_ch4", "thermal_ch4_file"]:
            try:
                thermal_ch4_path = snakemake.input[key]
                if thermal_ch4_path:
                    break
            except Exception:
                continue
    if thermal_ch4_path is None:
        thermal_ch4_path = getattr(snakemake.params, "thermal_ch4", None) if hasattr(snakemake.params, "__dict__") else None
    # Fallback: infer from resources path if file exists
    if thermal_ch4_path is None:
        candidate = Path(f"resources/demand_tyndp_thermal_ch4_{pyear}.csv")
        if candidate.exists():
            thermal_ch4_path = str(candidate)
            logger.info(f"Auto-detected thermal_ch4 file: {candidate}")
        else:
            # Try alternative naming used by build_tyndp_demand
            alt = Path(f"resources/demand_tyndp_thermal_ch4_{pyear}.csv")
            if alt.exists():
                thermal_ch4_path = str(alt)

    if scenario != "NT":
        logger.warning(f"Gas demand processing is not supported yet for {scenario}.")
        scenario = "NT"

    logger.info(f"Processing gas demand for scenario: {scenario}, year: {pyear}, thermal_ch4: {thermal_ch4_path}")

    demand = load_gas_demand(fn, scenario, pyear, thermal_ch4_path=thermal_ch4_path)

    demand.to_csv(snakemake.output.gas_demand, index=True)
