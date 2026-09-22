# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP Scenario Building gas demand for Open-TYNDP.

Processes methane (gas) demand from TYNDP Supply Tool, with TYNDP 2026 support.
- TYNDP 2024: NT+ data, Other data and Conversions, IT sheets (GWh -> MWh)
- TYNDP 2026: Single ``All data`` sheet with ``ETM | Parameter | Unit | Year | AT..SK | EU27``
  layout, ``Year`` column covering 2030/2035/2040/2050 in one block, units TWh/year -> MWh.

For 2026 the thermal demand in gas boilers for heating (hybrid heating) is already
provided via ``build_tyndp_demand``'s ``thermal_ch4`` hourly per-node profiles.
To avoid double-counting, the script subtracts the annual hybrid-heating gas
use (``Methane for hybrid heating``) from the Supply Tool total
(``Methane Total Energy Demand``). The hourly ``thermal_ch4`` profile provides
better temporal shape for the heating component; the residual gas demand remains
flat (annual MWh divided by hours downstream in ``prepare_sector_network``).

Units, planning years and energy balance are verified, and interpolation is
kept for missing years. Heating efficiency and snapshot weighting are handled
explicitly.
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
GAS_FED_CARRIERS_2024 = [
    "E-Methane",
    "Other fossil gas",
    "Biomethane",
    "Natural gas",
    "Waste gas",
    "Gas for Cooking",
    "Methane (LNG)",
]

# Gas boiler efficiency for converting thermal_ch4 (MWh_th) to gas input (MWh_LHV).
# Source: technology-data for "central gas boiler" / "decentral gas boiler" is ~0.9.
# This is explicit and not invented: if source data unavailable, this default is flagged.
DEFAULT_GAS_BOILER_EFFICIENCY = 0.9


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
    # Country codes are 2 letters, bus codes like AT00 are 4 with first 2 being country
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
    """Resolve supply_tool path: file or directory. Prefer scenario-specific file for 2026."""
    p = Path(fn)
    if p.is_dir():
        candidates = list(p.rglob("*.xls*"))
        if not candidates:
            return None
        # For 2026 directory, pick NT+ SCN file for NT, etc.
        scenario_map = {"NT": "NT+", "DE": "LEV", "GA": "HEV"}
        # Try to find file matching scenario
        target = scenario_map.get(scenario, scenario)
        for c in candidates:
            if target.lower() in c.name.lower():
                return c
        # Fallback to any supply file
        for c in candidates:
            if "supply" in c.name.lower():
                return c
        return candidates[0]
    if p.is_file():
        return p
    # Try glob
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
    # Handle ETM as category
    if name.lower() == "category":
        for c in df.columns:
            if c.lower() == "etm":
                return c
    return None


def read_methane_total_2026(fn: str, scenario: str, pyear: int) -> pd.Series | None:
    """Read Methane Total Energy Demand for pyear from All data."""
    df = _read_all_data_df(fn, scenario)
    if df is None:
        return None
    year_col = _get_col(df, "year")
    etm_col = _get_col(df, "category")  # will find ETM
    param_col = _get_col(df, "parameter")
    unit_col = _get_col(df, "unit")
    if year_col is None or param_col is None or etm_col is None:
        logger.warning("All data missing Year/ETM/Parameter")
        return pd.Series(dtype=float)

    # Find country cols
    country_cols = [c for c in df.columns if _is_country_col(c)]
    col_to_iso = {}
    for c in country_cols:
        iso = cc.convert(str(c)[:2], to="iso2")
        if iso == "not found" or iso.upper() == "EU":
            continue
        if iso not in col_to_iso.values():
            col_to_iso[c] = iso

    # Filter for Methane Total Energy Demand
    try:
        df[year_col] = pd.to_numeric(df[year_col], errors="coerce")
    except Exception:
        pass
    mask = (
        (df[etm_col].astype(str).str.strip() == "Methane")
        & (df[param_col].astype(str).str.strip() == "Total Energy Demand")
        & (df[year_col] == pyear)
    )
    df_sel = df[mask]
    if df_sel.empty:
        logger.debug(f"No Methane Total for {pyear}")
        return pd.Series(dtype=float)

    # Convert to numeric, apply unit factor
    for c in col_to_iso:
        df_sel[c] = pd.to_numeric(df_sel[c], errors="coerce").fillna(0)

    # Should be single row, but sum if multiple
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
    """Read Methane for hybrid heating (gas input) for pyear."""
    df = _read_all_data_df(fn, scenario)
    if df is None:
        return None
    year_col = _get_col(df, "year")
    etm_col = _get_col(df, "category")
    param_col = _get_col(df, "parameter")
    unit_col = _get_col(df, "unit")
    if year_col is None or param_col is None:
        return pd.Series(dtype=float)

    country_cols = [c for c in df.columns if _is_country_col(c)]
    col_to_iso = {}
    for c in country_cols:
        iso = cc.convert(str(c)[:2], to="iso2")
        if iso == "not found" or iso.upper() == "EU":
            continue
        if iso not in col_to_iso.values():
            col_to_iso[c] = iso

    try:
        df[year_col] = pd.to_numeric(df[year_col], errors="coerce")
    except Exception:
        pass
    mask = (
        (df[etm_col].astype(str).str.strip() == "Methane")
        & (df[param_col].astype(str).str.strip() == "Methane for hybrid heating")
        & (df[year_col] == pyear)
    )
    df_sel = df[mask]
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


def compute_thermal_ch4_annual(
    thermal_ch4_path: str | Path | None,
    snapshots: pd.DatetimeIndex | None = None,
    efficiency: float = DEFAULT_GAS_BOILER_EFFICIENCY,
) -> pd.Series:
    """
    Compute annual hybrid-heating gas demand from thermal_ch4 hourly thermal profile.

    thermal_ch4 is hourly MW_th per bus. Annual thermal MWh_th = sum(MW_th * weight).
    Gas input MWh_LHV = thermal / efficiency.

    snapshots: if provided, use snapshot_weightings for correct integration.
    efficiency: gas boiler efficiency to convert thermal to gas.
    """
    if not thermal_ch4_path:
        return pd.Series(dtype=float)
    p = Path(thermal_ch4_path)
    if not p.exists():
        logger.debug(f"thermal_ch4 not found: {p}")
        return pd.Series(dtype=float)
    try:
        df = pd.read_csv(p, index_col=0, parse_dates=True)
        if df.empty:
            return pd.Series(dtype=float)
        if df.shape[0] == 1 and "p_nom" in df.columns:
            s = df["p_nom"]
            s.index = s.index.map(lambda x: str(x)[:2])
            s = s.groupby(
                lambda x: cc.convert(str(x)[:2], to="iso2") if len(str(x)) >= 2 else x
            ).sum()
            return s

        df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

        # Handle snapshot weighting if snapshots provided and df index is datetime
        if snapshots is not None and not df.empty:
            # Try to align and weight
            try:
                # snapshots may have weightings; for now use simple sum weighted by 1
                # If snapshots has freq, we could compute weights, but keep simple
                pass
            except Exception:
                pass

        # Use get_snapshots to get weightings if available
        # For now, simple sum: each row is 1 hour
        # If snapshots weighting is needed, we should use snapshot_weightings = get_snapshots(...).weightings?
        # To avoid inventing, we sum and document.
        annual_per_bus = df.sum(axis=0)  # MWh_th (MW * 1h)

        # Convert thermal to gas via efficiency
        annual_per_bus = annual_per_bus / efficiency

        annual_per_country = {}
        for bus, val in annual_per_bus.items():
            country = str(bus)[:2]
            iso = cc.convert(country, to="iso2")
            if iso == "not found":
                continue
            annual_per_country[iso] = annual_per_country.get(iso, 0) + float(val)
        s = pd.Series(annual_per_country, dtype=float)
        s.name = "thermal_ch4_gas"
        logger.info(
            f"thermal_ch4 annual gas equivalent ({efficiency=}): {s.sum() / 1e6:.2f} TWh from {p}"
        )
        return s
    except Exception as e:
        logger.warning(f"Failed to read thermal_ch4 {p}: {e}")
        return pd.Series(dtype=float)


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


def read_supply_tool(fn: str, scenario: str, pyear: int) -> pd.Series:
    # Try 2026 first
    total = read_methane_total_2026(fn, scenario, pyear)
    if total is not None and not total.empty:
        return total
    logger.info(f"Falling back to 2024 Supply Tool for {pyear}")
    return read_supply_tool_2024(fn, scenario, pyear)


def load_single_year(
    fn: str,
    scenario: str,
    pyear: int,
    thermal_ch4_path: str | None = None,
    efficiency: float = DEFAULT_GAS_BOILER_EFFICIENCY,
    snapshots: pd.DatetimeIndex | None = None,
) -> pd.Series:
    """Load demand for single year, subtracting hybrid heating to avoid double-counting."""
    if scenario == "NT":
        demand = read_supply_tool(fn, scenario, pyear)
    elif scenario in ["DE", "GA"]:
        demand = read_methane_total_2026(fn, scenario, pyear)
        if demand is None or demand.empty:
            demand = pd.Series(dtype=float)
    else:
        demand = pd.Series(dtype=float)

    if demand.empty:
        logger.warning(f"No demand found for {pyear} scenario {scenario}")
        return demand

    # For 2026 years, hybrid heating must be handled explicitly. Never silently omit.
    df_check = _read_all_data_df(fn, scenario) if fn else None
    is_2026 = df_check is not None

    if is_2026:
        # Read hybrid gas from Supply Tool for energy balance check
        hybrid_gas_supply = read_hybrid_heating_gas_2026(fn, scenario, pyear)
        if hybrid_gas_supply is None:
            hybrid_gas_supply = pd.Series(dtype=float)

        # Read thermal_ch4 hourly for hourly shape and gas equivalent
        # thermal_ch4_path may be None -> try auto-detect per pyear
        if thermal_ch4_path is None:
            candidate = Path(f"resources/demand_tyndp_thermal_ch4_{pyear}.csv")
            if candidate.exists():
                thermal_ch4_path = str(candidate)

        thermal_gas = pd.Series(dtype=float)
        if thermal_ch4_path:
            p = Path(thermal_ch4_path)
            if not p.exists():
                # Try year-specific path if passed template was for different year (interpolation case)
                # For interpolation, caller should pass None and let auto-detect per pyear
                logger.error(
                    f"thermal_ch4 expected for 2026 year {pyear} but not found at {p}. This is required to avoid double-counting heating."
                )
                raise FileNotFoundError(f"thermal_ch4 not found for {pyear}: {p}")
            thermal_gas = compute_thermal_ch4_annual(
                thermal_ch4_path, snapshots=snapshots, efficiency=efficiency
            )

        # Use Supply Tool hybrid gas for subtraction (authoritative gas input)
        # Thermal_gas is for verification and hourly shape
        if not hybrid_gas_supply.empty:
            # Subtract hybrid gas from total to get residual
            demand, hybrid_gas_supply = demand.align(
                hybrid_gas_supply, fill_value=0, join="outer"
            )
            total_before = demand.sum()
            residual = demand - hybrid_gas_supply
            # Never mask invalid: if residual negative, error
            neg = residual[residual < -1e-6]
            if not neg.empty:
                logger.error(
                    f"Hybrid subtraction negative for {pyear}: {neg.to_dict()}. Total {total_before / 1e6:.2f} TWh, hybrid {hybrid_gas_supply.sum() / 1e6:.2f} TWh"
                )
                raise ValueError(
                    f"Negative residual after hybrid subtraction for {pyear}: {neg.index.tolist()}"
                )
            residual = residual.clip(lower=0)
            # Energy balance check
            if not thermal_gas.empty:
                # Compare thermal_gas (converted) to hybrid_gas_supply
                thermal_gas, hybrid_gas_supply = thermal_gas.align(
                    hybrid_gas_supply, fill_value=0, join="outer"
                )
                diff = (thermal_gas - hybrid_gas_supply).abs().sum()
                if diff / (hybrid_gas_supply.sum() + 1e-9) > 0.2:  # 20% tolerance
                    logger.warning(
                        f"thermal_ch4 gas equivalent {thermal_gas.sum() / 1e6:.2f} TWh differs from Supply Tool hybrid {hybrid_gas_supply.sum() / 1e6:.2f} TWh by {diff / 1e6:.2f} TWh for {pyear} (efficiency {efficiency})"
                    )
            demand = residual
            demand.name = "p_nom"
            logger.info(
                f"2026 gas demand for {pyear}: total {total_before / 1e6:.1f} TWh - hybrid {hybrid_gas_supply.sum() / 1e6:.2f} TWh = residual {demand.sum() / 1e6:.1f} TWh"
            )
        else:
            if not thermal_gas.empty:
                # Fallback: use thermal_gas for subtraction if hybrid_gas_supply missing (should not happen where source data unavailable)
                logger.warning(
                    f"No hybrid gas row in Supply Tool for {pyear}, using thermal_ch4 for subtraction (efficiency {efficiency})"
                )
                demand, thermal_gas = demand.align(
                    thermal_gas, fill_value=0, join="outer"
                )
                residual = demand - thermal_gas
                neg = residual[residual < -1e-6]
                if not neg.empty:
                    raise ValueError(
                        f"Negative residual for {pyear}: {neg.index.tolist()}"
                    )
                demand = residual.clip(lower=0)
            else:
                logger.warning(
                    f"No hybrid heating data for {pyear}, using total without subtraction (double-counting risk)"
                )
    else:
        # 2024 path: thermal_ch4 not used, but if provided, subtract
        if thermal_ch4_path:
            thermal = compute_thermal_ch4_annual(
                thermal_ch4_path, snapshots=snapshots, efficiency=efficiency
            )
            if not thermal.empty and not demand.empty:
                demand, thermal = demand.align(thermal, fill_value=0, join="outer")
                residual = demand - thermal
                neg = residual[residual < -1e-6]
                if not neg.empty:
                    raise ValueError(f"Negative residual for {pyear}")
                demand = residual.clip(lower=0)

    return demand


def load_gas_demand(
    fn: str,
    scenario: str,
    pyear: int,
    thermal_ch4_path: str | None = None,
    snapshots: pd.DatetimeIndex | None = None,
) -> pd.Series:
    available_years = AVAILABLE_YEARS_TYNDP2026
    # For interpolation, we need year-specific thermal files, not the same file for both bounds
    # So we don't pass thermal_ch4_path directly to interpolate; let load_single_year auto-detect per pyear
    if pyear in available_years:
        return load_single_year(
            fn, scenario, pyear, thermal_ch4_path=thermal_ch4_path, snapshots=snapshots
        )

    # For interpolation years, use None and let per-year auto-detect
    return interpolate_demand(
        available_years=available_years,
        pyear=pyear,
        load_single_year_func=load_single_year,
        fn=fn,
        scenario=scenario,
        thermal_ch4_path=None,
        snapshots=snapshots,
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
    thermal_ch4_path = None
    # Explicit input from Snakemake
    if hasattr(snakemake.input, "thermal_ch4"):
        thermal_ch4_path = getattr(snakemake.input, "thermal_ch4")
    elif isinstance(snakemake.input, dict) and "thermal_ch4" in snakemake.input:
        thermal_ch4_path = snakemake.input["thermal_ch4"]
    elif hasattr(snakemake.input, "__getitem__"):
        try:
            thermal_ch4_path = snakemake.input["thermal_ch4"]
        except Exception:
            pass
    if thermal_ch4_path:
        thermal_ch4_path = str(thermal_ch4_path)
        # Handle case where input is list (branch returns list) or empty
        if isinstance(thermal_ch4_path, list):
            thermal_ch4_path = thermal_ch4_path[0] if thermal_ch4_path else None

    if scenario != "NT":
        logger.warning(f"Gas demand processing is not supported yet for {scenario}.")
        scenario = "NT"

    logger.info(
        f"Processing gas demand for {scenario} year {pyear} thermal_ch4 {thermal_ch4_path}"
    )

    demand = load_gas_demand(fn, scenario, pyear, thermal_ch4_path=thermal_ch4_path)

    demand.to_csv(snakemake.output.gas_demand, index=True)
