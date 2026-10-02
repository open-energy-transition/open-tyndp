# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP Scenario Building gas demand for Open-TYNDP.

Processes methane (gas) demand from the TYNDP Supply Tool, with TYNDP 2026
support:

- TYNDP 2024: ``NT+ data``, ``Other data and Conversions`` and ``IT`` sheets
  of the single-file Supply Tool (GWh -> MWh). Unchanged legacy path.
- TYNDP 2026: single ``All data`` sheet with an ``ETM | Parameter | Unit |
  Year | AT..SK | EU27`` layout, ``Year`` column covering 2030/2035/2040/2050
  in one block, units TWh/year -> MWh.

For 2026 the gas demand for hybrid heating (``Methane for hybrid heating``)
is part of the Supply Tool total (``Methane Total Energy Demand``). Hourly
``thermal_ch4`` profiles built by ``build_tyndp_demand`` (upstream
``tyndp-2026``, PR #808) provide the temporal shape downstream; this script
outputs the residual annual gas demand (total minus hybrid gas) to avoid
double-counting, plus an hourly per-bus hybrid-heating series shaped by
``thermal_ch4`` and scaled to the Supply Tool hybrid totals. No
boiler-efficiency conversion is applied because both Supply Tool rows are
already gas input (MWh_LHV) and the thermal profiles are used as dimensionless
shape weights only.

Thermal-shape mapping (documented, not silently filled):

- Upstream ``read_demand_excel`` outputs long sheet-name columns
  (``Methane_Heat-<BUS>_HCH4-B``, UK already replaced by GB); bus codes are
  extracted with ``-<BUS>_`` and countries aggregate by first two letters.
- Upstream drops all-zero buses, so a country with hybrid gas but no thermal
  buses raises an explicit error instead of being silently flattened.
- Thermal files span 8760 hourly rows while snapshots may hold 8736h; the
  hourly output keeps the thermal index as-is and downstream reindexing
  (truncation/fill) belongs to ``prepare_sector_network``.
- Supply Tool ``All data`` has no GB column while thermal has GB/NI buses;
  multi-bus countries (ITx7, SEx4, DKx2) aggregate by first two letters.
  Remaining mapping gaps raise explicit errors; see #966.

Dependencies requiring maintainer direction:

- Scenario file mapping for the 2026 Supply-Tool directory (NT -> NT+) is
  assumed; DE -> LEV / GA -> HEV mapping is unconfirmed and logs a warning.
- Country coverage differences (e.g. GB/NI) are reported, not filled.
"""

import logging
import re
from pathlib import Path
from typing import Literal

import country_converter as coco
import pandas as pd

from scripts._helpers import configure_logging, interpolate_demand, set_scenario_config

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

    For a 2026 directory, prefer the scenario-specific file (NT -> NT+).
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


def _read_methane_rows(
    fn: str, scenario: str, planning_horizon: int, parameter: str
) -> pd.Series:
    """Read one Methane/parameter row selection for a year from All data."""
    df = _read_all_data_df(fn, scenario)
    if df is None:
        return pd.Series(dtype=float)
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
        & (df[param_col].astype(str).str.strip() == parameter)
        & (df[year_col] == planning_horizon)
    )
    df_sel = df[mask].copy()
    if df_sel.empty:
        logger.debug(f"No Methane '{parameter}' rows for {planning_horizon}")
        return pd.Series(dtype=float)

    for c in col_to_iso:
        df_sel[c] = pd.to_numeric(df_sel[c], errors="coerce").fillna(0)

    totals: dict[str, float] = {}
    for _, row in df_sel.iterrows():
        factor = (
            _unit_factor(str(row[unit_col]))
            if unit_col and unit_col in df_sel.columns
            else 1e6
        )
        for col, iso in col_to_iso.items():
            totals[iso] = totals.get(iso, 0.0) + float(row[col]) * factor

    s = pd.Series(totals, dtype=float)
    s = s[s != 0]
    logger.info(
        f"2026 Methane '{parameter}' for {planning_horizon}: {s.sum() / 1e6:.2f} TWh"
    )
    return s


def read_methane_total_2026(fn: str, scenario: str, planning_horizon: int) -> pd.Series:
    """Read Methane Total Energy Demand for a year from All data (MWh)."""
    s = _read_methane_rows(fn, scenario, planning_horizon, "Total Energy Demand")
    s.name = "p_nom"
    return s


def read_hybrid_heating_gas_2026(
    fn: str, scenario: str, planning_horizon: int
) -> pd.Series:
    """Read Methane for hybrid heating (gas input, MWh) for a year."""
    s = _read_methane_rows(fn, scenario, planning_horizon, "Methane for hybrid heating")
    s.name = "hybrid_gas"
    return s


def thermal_bus_of(column: str) -> str:
    """
    Extract the bus code from an upstream thermal demand column name.

    Upstream columns are long sheet names (``Methane_Heat-<BUS>_HCH4-B``);
    plain bus codes pass through unchanged.
    """
    m = re.search(r"-([A-Z0-9]+)_", str(column))
    return m.group(1) if m else str(column)


def read_thermal_hourly_csv(fn: str | Path) -> pd.DataFrame:
    """
    Read a ``build_tyndp_demand`` thermal_ch4 output CSV (hourly per bus).

    Returns shape weights indexed by datetime with one column per bus code.
    """
    df = pd.read_csv(fn, index_col=0, parse_dates=True)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    df.columns = [thermal_bus_of(c) for c in df.columns]
    # Merge duplicate bus columns (e.g. plain + corrected variants upstream).
    df = df.T.groupby(level=0).sum().T
    return df


def build_hybrid_hourly(
    thermal_df: pd.DataFrame, hybrid_annual: pd.Series
) -> pd.DataFrame:
    """
    Shape hourly per-bus hybrid heating gas demand from thermal profiles.

    Each country's buses are scaled so the country total equals the Supply
    Tool hybrid-heating annual for that country (shape from thermal, magnitude
    from gas input; no efficiency conversion invented). Buses in countries
    without Supply hybrid data are dropped with a warning. A country with
    hybrid gas but no thermal shape raises ``ValueError`` (data gap requiring
    maintainer direction, never silently flattened).
    """
    bus_annual = thermal_df.sum(axis=0)
    bus_country = {bus: str(bus)[:2] for bus in thermal_df.columns}
    country_thermal = {}
    for bus, val in bus_annual.items():
        country_thermal[bus_country[bus]] = country_thermal.get(
            bus_country[bus], 0.0
        ) + float(val)

    shaped = pd.DataFrame(0.0, index=thermal_df.index, columns=thermal_df.columns)
    for iso, target in hybrid_annual.items():
        iso = str(iso)
        buses = [b for b in thermal_df.columns if bus_country[b] == iso]
        country_shape = country_thermal.get(iso, 0.0)
        if not buses or country_shape <= 0:
            if float(target) > 0:
                raise ValueError(
                    f"Hybrid heating gas totals {float(target) / 1e6:.3f} TWh for "
                    f"{iso} but no thermal_ch4 shape buses (mapping gap, see #966)."
                )
            continue
        for bus in buses:
            shaped[bus] = thermal_df[bus] / country_shape * float(target)
    shaped = shaped.loc[:, (shaped.sum(axis=0) != 0)]
    # Energy-balance verification: shaped hourly sums equal hybrid annuals.
    for iso, target in hybrid_annual.items():
        iso = str(iso)
        buses = [b for b in shaped.columns if str(b)[:2] == iso]
        got = float(shaped[buses].sum().sum()) if buses else 0.0
        if abs(got - float(target)) > max(1.0, abs(float(target)) * 1e-9):
            raise ValueError(
                f"Hybrid shaping balance failed for {iso}: shaped {got:.1f} MWh "
                f"vs Supply hybrid {float(target):.1f} MWh."
            )
    shaped.columns.name = "Bus"
    shaped.index.name = shaped.index.name or "datetime"
    return shaped


def read_fed_data(
    fn: str, scenario: str, planning_horizon: int
) -> tuple[pd.Series, pd.Series]:
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
            skiprows=27 if planning_horizon == 2040 else 0,
            sheet_name="NT+ data",
        )

        # Set buses as column names
        demand_fed.columns = pd.Index(cc.convert(demand_fed.columns, to="iso2"))

        # Extract final heat demand
        demand_heat = demand_fed.loc["Heat"].mul(1e3)  # MWh

        # Extract final methane demand
        demand_fed = demand_fed.loc[GAS_FED_CARRIERS_2024].mul(1e3).sum()  # MWh

    except Exception as e:
        logger.warning(
            f"Failed to read final gas demand for scenario {scenario} and planning_horizon {planning_horizon}: "
            f"{type(e).__name__}: {e}"
        )
        demand_fed = pd.Series()
        demand_heat = pd.Series()

    return demand_fed, demand_heat


def read_heat_frame(
    fn: str, planning_horizon: int, type: Literal["distribution", "efficiency"]
) -> pd.DataFrame:
    """Read heat distribution and efficiency tables in 'Other data and Conversions' sheet of Supply Tool."""
    if type not in ["distribution", "efficiency"]:
        raise ValueError(
            f"Invalid type '{type}'. Must be 'distribution' or 'efficiency'."
        )
    if planning_horizon not in [2030, 2040, 2050]:
        raise ValueError(
            f"Invalid planning_horizon '{planning_horizon}'. Must be 2030, 2040, or 2050."
        )

    offset = ((planning_horizon - 2030) // 10) * 34
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
    df.columns = pd.Index(cc.convert(df.columns, to="iso2"))

    return df


def read_it_gas_prod(fn: str, planning_horizon: int) -> float:
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
        ).loc[planning_horizon, "For heat"]
        * 1e6
    )  # MWh


def read_heat_data(
    heat_fed: pd.Series, fn: str, scenario: str, planning_horizon: int
) -> pd.Series:
    """Read and process heat-related gas demand data from Supply Tool for a specific year."""
    try:
        shares = read_heat_frame(fn, planning_horizon, "distribution")
        efficiencies = read_heat_frame(fn, planning_horizon, "efficiency")

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
        demand.loc["IT"] = read_it_gas_prod(fn, planning_horizon)

    except Exception as e:
        logger.warning(
            f"Failed to read heat demand data for scenario {scenario} and planning_horizon {planning_horizon}: "
            f"{type(e).__name__}: {e}"
        )
        demand = pd.Series()

    return demand


def read_supply_tool_2024(fn: str, scenario: str, planning_horizon: int) -> pd.Series:
    """Read and process both final gas demand and heat demand data from Supply Tool for a specific year."""
    demand_fed, heat_fed = read_fed_data(fn, scenario, planning_horizon)
    demand_heat = read_heat_data(heat_fed, fn, scenario, planning_horizon)

    demand = pd.concat([demand_fed, demand_heat], axis=1).sum(axis=1)
    demand.name = "p_nom"

    return demand


def _is_2026_input(fn: str, scenario: str) -> bool:
    if not fn:
        return False
    return _read_all_data_df(fn, scenario) is not None


def demand_align(a: pd.Series, b: pd.Series):
    return a.align(b, fill_value=0, join="outer")


def load_single_year(fn: str, scenario: str, planning_horizon: int) -> pd.Series:
    """
    Load residual annual demand for a single year.

    For 2026: residual = Methane Total - Methane for hybrid heating.
    Raises ValueError on negative residual (invalid input, never masked).
    Raises ValueError when the total exists but the hybrid row is missing.
    For 2024: legacy FED + heat path (unchanged).
    """
    is_2026 = _is_2026_input(fn, scenario)

    if is_2026:
        if scenario in ["DE", "GA"]:
            logger.warning(
                "2026 DE/GA scenario file mapping (LEV/HEV) requires maintainer "
                f"confirmation; attempting scenario-specific file for {scenario}."
            )
        total = read_methane_total_2026(fn, scenario, planning_horizon)
        if total is None or total.empty:
            logger.warning(
                f"No 2026 Methane Total for {planning_horizon} scenario {scenario}"
            )
            return pd.Series(dtype=float, name="p_nom")

        hybrid = read_hybrid_heating_gas_2026(fn, scenario, planning_horizon)
        if hybrid is None or hybrid.empty:
            raise ValueError(
                f"2026 Methane Total exists for {planning_horizon} but "
                "'Methane for hybrid heating' row is missing; refusing to "
                "proceed without explicit heating data to avoid double-counting."
            )

        demand, hybrid = demand_align(total, hybrid)
        total_before = float(demand.sum())
        residual = demand - hybrid
        neg = residual[residual < -1e-6]
        if not neg.empty:
            raise ValueError(
                f"Negative residual after hybrid subtraction for {planning_horizon}: "
                f"{neg.index.tolist()}. Total {total_before / 1e6:.2f} TWh, "
                f"hybrid {float(hybrid.sum()) / 1e6:.2f} TWh."
            )
        residual = residual.clip(lower=0)
        residual.name = "p_nom"
        # Energy balance verification (exact, no tolerance threshold)
        balance = float(residual.sum() + hybrid.sum())
        logger.info(
            f"2026 gas demand for {planning_horizon}: total {total_before / 1e6:.1f} TWh - "
            f"hybrid {float(hybrid.sum()) / 1e6:.2f} TWh = "
            f"residual {float(residual.sum() / 1e6):.1f} TWh "
            f"(balance {balance / 1e6:.1f} TWh)"
        )
        return residual

    # 2024 legacy path (unchanged, preserves existing workflow)
    if scenario == "NT":
        return read_supply_tool_2024(fn, scenario, planning_horizon)
    logger.warning(f"Gas demand processing is not supported yet for {scenario}.")
    return pd.Series(dtype=float, name="p_nom")


def load_gas_demand(fn: str, scenario: str, planning_horizon: int) -> pd.Series:
    """Load gas demand, selecting 2026 vs 2024 available years by input type."""
    is_2026 = _is_2026_input(fn, scenario)
    available_years = (
        AVAILABLE_YEARS_TYNDP2026 if is_2026 else AVAILABLE_YEARS_TYNDP2024
    )

    if planning_horizon in available_years:
        return load_single_year(fn, scenario, planning_horizon)

    return interpolate_demand(
        available_years=available_years,
        planning_horizon=planning_horizon,
        load_single_year_func=load_single_year,
        fn=fn,
        scenario=scenario,
    )


def load_hybrid_hourly(
    fn: str, scenario: str, planning_horizon: int, thermal_fn: str | Path | None
) -> pd.DataFrame:
    """
    Build hourly per-bus hybrid heating gas demand for a 2026 year.

    Returns an empty DataFrame when the input is not 2026 or no thermal
    file is provided (2024 behaviour: no hybrid series).
    """
    if not thermal_fn or not _is_2026_input(fn, scenario):
        return pd.DataFrame()
    if isinstance(thermal_fn, (list, tuple)):
        thermal_fn = thermal_fn[0] if thermal_fn else None
    if thermal_fn is None:
        return pd.DataFrame()
    hybrid = read_hybrid_heating_gas_2026(fn, scenario, planning_horizon)
    if hybrid is None or hybrid.empty:
        raise ValueError(
            f"2026 Methane Total exists for {planning_horizon} but "
            "'Methane for hybrid heating' row is missing; refusing to "
            "proceed without explicit heating data to avoid double-counting."
        )
    thermal_df = read_thermal_hourly_csv(thermal_fn)
    if thermal_df.empty:
        raise ValueError(
            f"Thermal_ch4 hourly input {thermal_fn} holds no demand; cannot "
            "shape hybrid heating (mapping gap, see #966)."
        )
    return build_hybrid_hourly(thermal_df, hybrid)


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

    # Parameters
    scenario = snakemake.params["scenario"]
    fn = snakemake.input.supply_tool
    planning_horizon = int(snakemake.wildcards.planning_horizons)
    thermal_fn = snakemake.input.get("thermal_ch4", None)

    # Load residual annual demand with interpolation
    logger.info(f"Processing gas demand for scenario: {scenario}")
    demand = load_gas_demand(fn, scenario, planning_horizon)

    # Hourly per-bus hybrid heating shaped by thermal_ch4 (2026 only)
    try:
        hybrid_hourly = load_hybrid_hourly(fn, scenario, planning_horizon, thermal_fn)
    except ValueError:
        raise

    # Export to CSV
    demand.to_csv(snakemake.output.gas_demand, index=True)
    if hybrid_hourly.empty:
        pd.DataFrame(columns=["Bus"]).to_csv(snakemake.output.gas_hybrid, index=True)
        logger.warning(
            "No hybrid hourly series written (2024 input or no thermal_ch4); "
            "wrote empty hybrid file."
        )
    else:
        hybrid_hourly.to_csv(snakemake.output.gas_hybrid, index=True)
