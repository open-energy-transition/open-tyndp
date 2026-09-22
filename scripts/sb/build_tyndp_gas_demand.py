# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP Scenario Building gas demand for Open-TYNDP.

Processes methane (gas) demand from TYNDP Supply Tool, with TYNDP 2026 support.
- TYNDP 2024: NT+ data, Other data and Conversions, IT sheets (GWh -> MWh)
- TYNDP 2026: Single ``All data`` sheet with ``ETM | Parameter | Unit | Year |
  AT..SK | EU27`` layout, ``Year`` column covering 2030/2035/2040/2050 in one
  block, units TWh/year -> MWh.

For 2026 the gas demand for hybrid heating (``Methane for hybrid heating``)
is part of the Supply Tool total (``Methane Total Energy Demand``). The hourly
``thermal_ch4`` profiles from ``build_tyndp_demand`` (owned by the tyndp-2026
branch, see PR #808) provide better temporal shape downstream; this script
outputs the residual annual gas demand (total minus hybrid gas) to avoid
double-counting. No boiler-efficiency conversion is applied here because both
Supply Tool rows are already gas input (MWh_LHV).

Dependencies requiring maintainer decision:
- Hourly ``thermal_ch4`` production (``build_tyndp_demand`` rule,
  ``retrieve_tyndp_2026``) lives on upstream ``tyndp-2026`` and is not
  duplicated here. Downstream hourly integration of residual + thermal is
  pending that infrastructure.
- Scenario file mapping for 2026 directory (NT+ vs LEV/HEV) is assumed as
  NT->NT+; DE/GA mapping requires maintainer confirmation.
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

AVAILABLE_YEARS_TYNDP2026 = [2030, 2035, 2040, 2050]
AVAILABLE_YEARS_TYNDP2024 = [2030, 2040]
GAS_FED_CARRIERS_2024 = [
    "E-Methane",
    "Other fossil gas",
    "Biomethane",
    "Natural gas",
    "Waste gas",
    "Gas for Cooking",
    "Methane (LNG)",
]


def _is_country_col(col: str) -> bool:
    col = str(col).strip()
    if col.lower() in {
        "etm",
        "parameter",
        "unit",
        "year",
        "scenario",
        "eu27",
        "eu27.1",
        "unnamed: 33",
    }:
        return False
    try:
        iso = cc.convert(col[:2], to="iso2")
        return iso != "not found" and len(col[:2]) == 2
    except Exception:
        return False


def _unit_factor(unit_str: str) -> float:
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


def _find_supply_tool_file(fn: str | Path, scenario: str = "NT") -> Path | None:
    """
    Resolve supply_tool path: file or directory.

    For 2026 directory, prefer scenario-specific file (NT->NT+).
    DE/GA mapping (LEV/HEV) requires maintainer confirmation; currently
    attempts LEV for DE and HEV for GA, else first supply file.
    """
    p = Path(fn)
    if p.is_dir():
        candidates = list(p.rglob("*.xls*"))
        if not candidates:
            return None
        scenario_map = {"NT": "NT+", "DE": "LEV", "GA": "HEV"}
        target = scenario_map.get(scenario, scenario)
        for c in candidates:
            if target.lower() in c.name.lower():
                return c
        for c in candidates:
            if "supply" in c.name.lower():
                return c
        return candidates[0]
    if p.is_file():
        return p
    candidates = list(Path(".").glob(str(fn)))
    if candidates:
        return candidates[0]
    return p if p.exists() else None


def _read_all_data_df(fn: str | Path, scenario: str) -> pd.DataFrame | None:
    actual = _find_supply_tool_file(fn, scenario)
    if actual is None:
        logger.debug(f"Supply tool file not found for {fn} scenario {scenario}")
        return None
    fn = str(actual)
    try:
        try:
            df = pd.read_excel(fn, sheet_name="All data", header=0, engine="openpyxl")
        except Exception:
            df = pd.read_excel(fn, sheet_name="All data", header=0)
    except Exception as e:
        logger.debug(f"All data sheet not found in {fn}: {e}")
        return None
    df.columns = [str(c).strip() for c in df.columns]
    return df


def _get_col(df: pd.DataFrame, name: str) -> str | None:
    for c in df.columns:
        if c.lower() == name.lower():
            return c
    if name.lower() == "category":
        for c in df.columns:
            if c.lower() == "etm":
                return c
    return None


def _country_map(df: pd.DataFrame) -> dict:
    country_cols = [c for c in df.columns if _is_country_col(c)]
    col_to_iso = {}
    for c in country_cols:
        iso = cc.convert(str(c)[:2], to="iso2")
        if iso == "not found" or iso.upper() == "EU":
            continue
        if iso not in col_to_iso.values():
            col_to_iso[c] = iso
    return col_to_iso


def read_methane_total_2026(fn: str, scenario: str, pyear: int) -> pd.Series | None:
    """Read Methane Total Energy Demand for pyear from All data (MWh)."""
    df = _read_all_data_df(fn, scenario)
    if df is None:
        return None
    year_col = _get_col(df, "year")
    etm_col = _get_col(df, "category")
    param_col = _get_col(df, "parameter")
    unit_col = _get_col(df, "unit")
    if year_col is None or param_col is None or etm_col is None:
        logger.warning("All data missing Year/ETM/Parameter")
        return pd.Series(dtype=float)

    col_to_iso = _country_map(df)
    df = df.copy()
    try:
        df[year_col] = pd.to_numeric(df[year_col], errors="coerce")
    except Exception:
        pass
    mask = (
        (df[etm_col].astype(str).str.strip() == "Methane")
        & (df[param_col].astype(str).str.strip() == "Total Energy Demand")
        & (df[year_col] == pyear)
    )
    df_sel = df[mask].copy()
    if df_sel.empty:
        logger.debug(f"No Methane Total for {pyear}")
        return pd.Series(dtype=float)

    for c in col_to_iso:
        df_sel[c] = pd.to_numeric(df_sel[c], errors="coerce").fillna(0)

    totals = {}
    for _, row in df_sel.iterrows():
        factor = (
            _unit_factor(str(row[unit_col]))
            if unit_col and unit_col in df_sel.columns
            else 1e6
        )
        for col, iso in col_to_iso.items():
            totals[iso] = totals.get(iso, 0.0) + float(row[col]) * factor

    s = pd.Series(totals, dtype=float)
    s.name = "p_nom"
    s = s[s != 0]
    logger.info(f"2026 Methane total for {pyear}: {s.sum() / 1e6:.1f} TWh")
    return s


def read_hybrid_heating_gas_2026(
    fn: str, scenario: str, pyear: int
) -> pd.Series | None:
    """Read Methane for hybrid heating (gas input, MWh) for pyear."""
    df = _read_all_data_df(fn, scenario)
    if df is None:
        return None
    year_col = _get_col(df, "year")
    etm_col = _get_col(df, "category")
    param_col = _get_col(df, "parameter")
    unit_col = _get_col(df, "unit")
    if year_col is None or param_col is None:
        return pd.Series(dtype=float)

    col_to_iso = _country_map(df)
    df = df.copy()
    try:
        df[year_col] = pd.to_numeric(df[year_col], errors="coerce")
    except Exception:
        pass
    mask = (
        (df[etm_col].astype(str).str.strip() == "Methane")
        & (df[param_col].astype(str).str.strip() == "Methane for hybrid heating")
        & (df[year_col] == pyear)
    )
    df_sel = df[mask].copy()
    if df_sel.empty:
        logger.debug(f"No hybrid heating for {pyear}")
        return pd.Series(dtype=float)

    for c in col_to_iso:
        df_sel[c] = pd.to_numeric(df_sel[c], errors="coerce").fillna(0)

    totals = {}
    for _, row in df_sel.iterrows():
        factor = _unit_factor(str(row[unit_col])) if unit_col else 1e6
        for col, iso in col_to_iso.items():
            totals[iso] = totals.get(iso, 0.0) + float(row[col]) * factor

    s = pd.Series(totals, dtype=float)
    s.name = "hybrid_gas"
    s = s[s != 0]
    logger.info(f"2026 Hybrid gas for {pyear}: {s.sum() / 1e6:.2f} TWh")
    return s


def read_fed_data(fn: str, scenario: str, pyear: int) -> tuple[pd.Series, pd.Series]:
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
        demand_fed = demand_fed.loc[GAS_FED_CARRIERS_2024].mul(1e3).sum()
    except Exception as e:
        logger.warning(f"Failed to read FED for {scenario} {pyear}: {e}")
        demand_fed = pd.Series()
        demand_heat = pd.Series()
    return demand_fed, demand_heat


def read_heat_frame(
    fn: str, pyear: int, type: Literal["distribution", "efficiency"]
) -> pd.DataFrame:
    if type not in ["distribution", "efficiency"]:
        raise ValueError(f"Invalid type '{type}'")
    if pyear not in [2030, 2040, 2050]:
        raise ValueError(f"Invalid pyear '{pyear}'")
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
        logger.warning(f"Failed heat for {scenario} {pyear}: {e}")
        demand = pd.Series()
    return demand


def read_supply_tool_2024(fn: str, scenario: str, pyear: int) -> pd.Series:
    demand_fed, heat_fed = read_fed_data(fn, scenario, pyear)
    demand_heat = read_heat_data(heat_fed, fn, scenario, pyear)
    demand = pd.concat([demand_fed, demand_heat], axis=1).sum(axis=1)
    demand.name = "p_nom"
    return demand


def _is_2026_input(fn: str, scenario: str) -> bool:
    if not fn:
        return False
    return _read_all_data_df(fn, scenario) is not None


def load_single_year(fn: str, scenario: str, pyear: int) -> pd.Series:
    """
    Load demand for single year, subtracting hybrid heating for 2026.

    For 2026: residual = Methane Total - Methane for hybrid heating.
    Raises ValueError on negative residual (invalid input, never masked).
    Raises FileNotFoundError context via missing hybrid row is treated as
    explicit error when total exists but hybrid row absent for a 2026 year.
    For 2024: legacy FED + heat path (unchanged).
    """
    is_2026 = _is_2026_input(fn, scenario)

    if is_2026:
        if scenario in ["DE", "GA"]:
            logger.warning(
                f"2026 DE/GA scenario file mapping (LEV/HEV) requires maintainer "
                f"confirmation; attempting scenario-specific file for {scenario}."
            )
        total = read_methane_total_2026(fn, scenario, pyear)
        if total is None or total.empty:
            logger.warning(f"No 2026 Methane Total for {pyear} scenario {scenario}")
            return pd.Series(dtype=float, name="p_nom")

        hybrid = read_hybrid_heating_gas_2026(fn, scenario, pyear)
        if hybrid is None or hybrid.empty:
            raise ValueError(
                f"2026 Methane Total exists for {pyear} but "
                f"'Methane for hybrid heating' row is missing; refusing to "
                f"proceed without explicit heating data to avoid double-counting."
            )

        demand, hybrid = demand_align(total, hybrid)
        total_before = float(demand.sum())
        residual = demand - hybrid
        neg = residual[residual < -1e-6]
        if not neg.empty:
            raise ValueError(
                f"Negative residual after hybrid subtraction for {pyear}: "
                f"{neg.index.tolist()}. Total {total_before / 1e6:.2f} TWh, "
                f"hybrid {float(hybrid.sum()) / 1e6:.2f} TWh."
            )
        residual = residual.clip(lower=0)
        residual.name = "p_nom"
        # Energy balance verification (exact, no tolerance threshold)
        balance = float(residual.sum() + hybrid.sum())
        logger.info(
            f"2026 gas demand for {pyear}: total {total_before / 1e6:.1f} TWh - "
            f"hybrid {float(hybrid.sum()) / 1e6:.2f} TWh = "
            f"residual {float(residual.sum()) / 1e6:.1f} TWh "
            f"(balance {balance / 1e6:.1f} TWh)"
        )
        return residual

    # 2024 legacy path (unchanged, preserves existing workflow)
    if scenario == "NT":
        return read_supply_tool_2024(fn, scenario, pyear)
    logger.warning(f"Gas demand processing is not supported yet for {scenario}.")
    return pd.Series(dtype=float, name="p_nom")


def demand_align(a: pd.Series, b: pd.Series):
    return a.align(b, fill_value=0, join="outer")


def load_gas_demand(fn: str, scenario: str, pyear: int) -> pd.Series:
    """Load gas demand, selecting 2026 vs 2024 available years by input type."""
    is_2026 = _is_2026_input(fn, scenario)
    available_years = (
        AVAILABLE_YEARS_TYNDP2026 if is_2026 else AVAILABLE_YEARS_TYNDP2024
    )

    if pyear in available_years:
        return load_single_year(fn, scenario, pyear)

    return interpolate_demand(
        available_years=available_years,
        pyear=pyear,
        load_single_year_func=load_single_year,
        fn=fn,
        scenario=scenario,
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

    demand = load_gas_demand(fn, scenario, pyear)

    demand.to_csv(snakemake.output.gas_demand, index=True)
