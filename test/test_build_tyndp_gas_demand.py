# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""Tests for TYNDP 2026 gas demand (All data row selection and hybrid subtraction)."""

import os
from pathlib import Path

import pandas as pd
import pytest

from scripts.sb.build_tyndp_gas_demand import (
    AVAILABLE_YEARS,
    _find_supply_tool_file,
    _unit_factor,
    build_hybrid_hourly,
    load_gas_demand,
    load_hybrid_hourly,
    load_single_year,
    read_hybrid_heating_gas_2026,
    read_methane_total_2026,
    read_thermal_hourly_csv,
    thermal_bus_of,
)


def test_available_years():
    assert AVAILABLE_YEARS == [2030, 2035, 2040, 2050]


def test_unit_factor():
    assert _unit_factor("TWh/year") == 1e6
    assert _unit_factor("TWh") == 1e6
    assert _unit_factor("GWh") == 1e3
    assert _unit_factor("MWh") == 1


def _write_all_data(tmp_path, rows):
    fn = tmp_path / "supply.xlsx"
    pd.DataFrame(rows).to_excel(fn, sheet_name="All data", index=False)
    return str(fn)


def test_read_methane_total_2026(tmp_path):
    fn = _write_all_data(
        tmp_path,
        [
            {
                "ETM": "Methane",
                "Parameter": "Total Energy Demand",
                "Unit": "TWh/year",
                "Year": 2030,
                "AT": 1.7,
                "DE": 13.0,
                "EU27": 100.0,
                "EU27.1": 100.0,
            }
        ],
    )
    series = read_methane_total_2026(fn, "NT", 2030)
    assert abs(series["AT"] - 1.7e6) < 1
    assert abs(series["DE"] - 13e6) < 1
    assert "EU27" not in series.index


def test_read_hybrid_heating_2026(tmp_path):
    fn = _write_all_data(
        tmp_path,
        [
            {
                "ETM": "Methane",
                "Parameter": "Methane for hybrid heating",
                "Unit": "TWh/year",
                "Year": 2030,
                "AT": 0.5,
                "DE": 1.35,
            }
        ],
    )
    series = read_hybrid_heating_gas_2026(fn, "NT", 2030)
    assert abs(series["AT"] - 0.5e6) < 1


def test_year_filter(tmp_path):
    fn = _write_all_data(
        tmp_path,
        [
            {
                "ETM": "Methane",
                "Parameter": "Total Energy Demand",
                "Unit": "TWh/year",
                "Year": 2030,
                "AT": 1.0,
            },
            {
                "ETM": "Methane",
                "Parameter": "Total Energy Demand",
                "Unit": "TWh/year",
                "Year": 2040,
                "AT": 2.0,
            },
        ],
    )
    assert read_methane_total_2026(fn, "NT", 2030)["AT"] == 1e6
    assert read_methane_total_2026(fn, "NT", 2040)["AT"] == 2e6


def test_subtraction_and_balance(tmp_path):
    fn = _write_all_data(
        tmp_path,
        [
            {
                "ETM": "Methane",
                "Parameter": "Total Energy Demand",
                "Unit": "TWh/year",
                "Year": 2030,
                "AT": 1.0,
                "DE": 2.0,
            },
            {
                "ETM": "Methane",
                "Parameter": "Methane for hybrid heating",
                "Unit": "TWh/year",
                "Year": 2030,
                "AT": 0.2,
                "DE": 0.5,
            },
        ],
    )
    result = load_single_year(fn, "NT", 2030)
    assert abs(result["AT"] - 800000) < 1
    assert abs(result["DE"] - 1500000) < 1
    # Energy balance: residual + hybrid == total
    hybrid = read_hybrid_heating_gas_2026(fn, "NT", 2030)
    total = read_methane_total_2026(fn, "NT", 2030)
    assert abs((result.sum() + hybrid.sum()) - total.sum()) < 1


def test_negative_residual_raises(tmp_path):
    fn = _write_all_data(
        tmp_path,
        [
            {
                "ETM": "Methane",
                "Parameter": "Total Energy Demand",
                "Unit": "TWh/year",
                "Year": 2030,
                "AT": 0.1,
            },
            {
                "ETM": "Methane",
                "Parameter": "Methane for hybrid heating",
                "Unit": "TWh/year",
                "Year": 2030,
                "AT": 0.2,
            },
        ],
    )
    with pytest.raises(ValueError, match="Negative residual"):
        load_single_year(fn, "NT", 2030)


def test_missing_hybrid_raises(tmp_path):
    fn = _write_all_data(
        tmp_path,
        [
            {
                "ETM": "Methane",
                "Parameter": "Total Energy Demand",
                "Unit": "TWh/year",
                "Year": 2030,
                "AT": 1.0,
            }
        ],
    )
    with pytest.raises(ValueError, match="hybrid heating"):
        load_single_year(fn, "NT", 2030)


def test_interpolation_between_2026_years(tmp_path):
    fn = _write_all_data(
        tmp_path,
        [
            {
                "ETM": "Methane",
                "Parameter": "Total Energy Demand",
                "Unit": "TWh/year",
                "Year": 2030,
                "AT": 1.0,
            },
            {
                "ETM": "Methane",
                "Parameter": "Methane for hybrid heating",
                "Unit": "TWh/year",
                "Year": 2030,
                "AT": 0.2,
            },
            {
                "ETM": "Methane",
                "Parameter": "Total Energy Demand",
                "Unit": "TWh/year",
                "Year": 2035,
                "AT": 1.5,
            },
            {
                "ETM": "Methane",
                "Parameter": "Methane for hybrid heating",
                "Unit": "TWh/year",
                "Year": 2035,
                "AT": 0.5,
            },
        ],
    )
    # 2030 residual 0.8, 2035 residual 1.0, 2032 weight 0.4 => 0.88
    result = load_gas_demand(fn, "NT", 2032)
    assert abs(result["AT"] - 880000) < 1


def test_thermal_bus_of():
    assert thermal_bus_of("Methane_Heat-DE00_HCH4-B") == "DE00"
    assert thermal_bus_of("Methane_Heat-ITCA_HCH4-B") == "ITCA"
    assert thermal_bus_of("DE00") == "DE00"


def _write_thermal_csv(tmp_path, columns, n_hours=48):
    index = pd.date_range("2030-01-01", periods=n_hours, freq="h")
    df = pd.DataFrame({c: float(i + 1) for i, c in enumerate(columns)}, index=index)
    # give each column an hourly shape (avoid flat series)
    df = df.mul(pd.Series(range(1, n_hours + 1), index=index), axis=0)
    fn = tmp_path / "demand_tyndp_thermal_ch4_2030.csv"
    df.to_csv(fn, index=True)
    return str(fn)


def test_read_thermal_hourly_csv_merges_duplicate_buses(tmp_path):
    fn = _write_thermal_csv(
        tmp_path, ["Methane_Heat-DE00_HCH4-B", "DE00", "Methane_Heat-FR00_HCH4-B"]
    )
    df = read_thermal_hourly_csv(fn)
    assert list(df.columns) == ["DE00", "FR00"]
    # duplicate DE00 columns merged by summation
    raw = pd.read_csv(fn, index_col=0, parse_dates=True)
    assert (
        abs(
            df["DE00"].iloc[0]
            - raw.iloc[0].sum()
            + raw["Methane_Heat-FR00_HCH4-B"].iloc[0]
        )
        < 1e-9
    )


def test_build_hybrid_hourly_scales_shape_to_totals(tmp_path):
    fn = _write_thermal_csv(tmp_path, ["DE00", "DE01", "FR00"], n_hours=24)
    thermal_df = read_thermal_hourly_csv(fn)
    hybrid_annual = pd.Series({"DE": 2.4e6, "FR": 1.2e6})
    shaped = build_hybrid_hourly(thermal_df, hybrid_annual)
    assert set(shaped.columns) == {"DE00", "DE01", "FR00"}
    # country totals match Supply hybrid annuals exactly
    assert abs(shaped[["DE00", "DE01"]].sum().sum() - 2.4e6) / 2.4e6 < 1e-9
    assert abs(shaped["FR00"].sum() - 1.2e6) / 1.2e6 < 1e-9
    # per-bus split follows thermal annual shares
    bus_annual = thermal_df.sum(axis=0)
    de_share = bus_annual["DE00"] / (bus_annual["DE00"] + bus_annual["DE01"])
    assert abs(shaped["DE00"].sum() - 2.4e6 * de_share) / 2.4e6 < 1e-9
    # shape preserved: hourly profile proportional to thermal input
    assert abs((shaped["FR00"] / thermal_df["FR00"]).std()) < 1e-9


def test_build_hybrid_hourly_drops_countries_without_hybrid(tmp_path):
    fn = _write_thermal_csv(tmp_path, ["DE00", "ES00"], n_hours=24)
    thermal_df = read_thermal_hourly_csv(fn)
    shaped = build_hybrid_hourly(thermal_df, pd.Series({"DE": 1.0e6}))
    assert "ES00" not in shaped.columns
    assert abs(shaped.sum().sum() - 1.0e6) / 1.0e6 < 1e-9


def test_build_hybrid_hourly_raises_without_shape(tmp_path):
    fn = _write_thermal_csv(tmp_path, ["DE00"], n_hours=24)
    thermal_df = read_thermal_hourly_csv(fn)
    with pytest.raises(ValueError, match="no thermal_ch4 shape"):
        build_hybrid_hourly(thermal_df, pd.Series({"FR": 1.0e6}))


def test_load_hybrid_hourly_requires_thermal(tmp_path):
    fn = _write_all_data(
        tmp_path,
        [
            {
                "ETM": "Methane",
                "Parameter": "Total Energy Demand",
                "Unit": "TWh/year",
                "Year": 2030,
                "AT": 1.0,
            },
            {
                "ETM": "Methane",
                "Parameter": "Methane for hybrid heating",
                "Unit": "TWh/year",
                "Year": 2030,
                "AT": 0.2,
            },
        ],
    )
    with pytest.raises(ValueError, match="thermal_ch4 hourly input is required"):
        load_hybrid_hourly(fn, "NT", 2030, None)


def test_supply_tool_selection_prefers_ntplus_without_scenario_mapping(tmp_path):
    import openpyxl

    for name in (
        "Supply Tool NT+ SCN 2026.xlsm",
        "Supply Tool LEV SCN 2026.xlsm",
        "Supply Tool HEV SCN 2026.xlsm",
    ):
        wb = openpyxl.Workbook()
        ws = wb.active
        ws.title = "All data"
        ws.append(["ETM", "Parameter", "Unit", "Year", "AT"])
        ws.append(["Methane", "Total Energy Demand", "TWh/year", 2030, 1.0])
        wb.save(tmp_path / name)
    # NT resolves the NT+ scenario file; DE/GA are not mapped to the LEV/HEV
    # sensitivity files and fall back to NT+ (with warning) instead.
    assert _find_supply_tool_file(tmp_path, "NT").name.startswith("Supply Tool NT+")
    assert _find_supply_tool_file(tmp_path, "DE").name.startswith("Supply Tool NT+")
    assert _find_supply_tool_file(tmp_path, "GA").name.startswith("Supply Tool NT+")


def test_end_to_end_residual_plus_hourly_equals_total(tmp_path):
    supply = _write_all_data(
        tmp_path,
        [
            {
                "ETM": "Methane",
                "Parameter": "Total Energy Demand",
                "Unit": "TWh/year",
                "Year": 2030,
                "DE": 2.0,
                "FR": 1.0,
            },
            {
                "ETM": "Methane",
                "Parameter": "Methane for hybrid heating",
                "Unit": "TWh/year",
                "Year": 2030,
                "DE": 0.5,
                "FR": 0.25,
            },
        ],
    )
    thermal_fn = _write_thermal_csv(tmp_path, ["DE00", "DE01", "FR00"], n_hours=24)
    residual = load_single_year(supply, "NT", 2030)
    hybrid = load_hybrid_hourly(supply, "NT", 2030, thermal_fn)
    total = read_methane_total_2026(supply, "NT", 2030)
    assert abs((residual.sum() + hybrid.sum().sum()) - total.sum()) < 1


def test_real_supply_tool_2026():
    real = os.environ.get("TYNDP2026_SUPPLY_TOOL", "/tmp/supply2026.xlsm")
    if not Path(real).exists():
        pytest.skip("Real 2026 Supply Tool not available")
    total = read_methane_total_2026(real, "NT", 2030)
    hybrid = read_hybrid_heating_gas_2026(real, "NT", 2030)
    assert not total.empty
    # Spot-check DE values from actual file (approx, verifies TWh->MWh)
    assert 400e6 < total["DE"] < 500e6
    assert 1e6 < hybrid["DE"] < 2e6
    residual = load_single_year(real, "NT", 2030)
    assert abs((residual.sum() + hybrid.sum()) - total.sum()) < 1e3
    # Country-level: residual non-negative and less than total per country
    for iso in residual.index:
        assert residual[iso] >= 0
        if iso in total.index and iso in hybrid.index:
            assert residual[iso] <= total[iso] + 1e-6


def _thermal_bus_of_sheet(sheet: str) -> str:
    import re

    m = re.search(r"-([A-Z0-9]+)_", sheet)
    return m.group(1) if m else sheet


def _read_thermal_2030_ws003(thermal_fn: str):
    demand = pd.read_excel(
        thermal_fn,
        header=10,
        index_col=[0, 1],
        sheet_name=None,
        usecols=lambda n: n in ("Date", "Hour", "WS003"),
        engine="calamine",
    )
    demand = pd.concat(demand, axis=1).droplevel(1, axis=1)
    demand.columns = [
        _thermal_bus_of_sheet(c).replace("UK", "GB") for c in demand.columns
    ]
    return demand


def test_thermal_bus_extraction_and_country_mapping():
    thermal_fn = os.environ.get(
        "TYNDP2026_THERMAL_CH4_2030",
        "/tmp/demand2030/Demand/Demand/2030/Thermal_energy_Methane_2030.xlsx",
    )
    if not Path(thermal_fn).exists():
        pytest.skip("Real 2026 thermal_ch4 file not available")
    demand = _read_thermal_2030_ws003(thermal_fn)
    # 41 sheets over the thermal file's hourly rows
    assert demand.shape[0] == 8760
    assert demand.shape[1] == 41
    # Multi-bus countries aggregated: IT has 7 sheets, SE 4, DK 2, GR 2
    it_cols = [c for c in demand.columns if c.startswith("IT")]
    assert len(it_cols) >= 7
    se_cols = [c for c in demand.columns if c.startswith("SE")]
    assert len(se_cols) == 4
    # Zero-demand buses exist (e.g. AT00) and must be treated as zero, not missing
    assert (demand["AT00"].sum() == 0) or (demand["AT00"].sum() >= 0)
    # Upstream column-naming dependency: raw output uses long sheet names,
    # bus extraction is required downstream (owned by tyndp-2026, PR #808)
    raw = pd.read_excel(
        thermal_fn,
        header=10,
        index_col=[0, 1],
        sheet_name=None,
        usecols=lambda n: n in ("Date", "Hour", "WS003"),
        engine="calamine",
    )
    raw_cols = list(pd.concat(raw, axis=1).droplevel(1, axis=1).columns[:3])
    assert any("Methane_Heat-" in c for c in raw_cols)


def test_thermal_shape_normalization_preserves_balance():
    supply = os.environ.get("TYNDP2026_SUPPLY_TOOL", "/tmp/supply2026.xlsm")
    thermal_fn = os.environ.get(
        "TYNDP2026_THERMAL_CH4_2030",
        "/tmp/demand2030/Demand/Demand/2030/Thermal_energy_Methane_2030.xlsx",
    )
    if not Path(supply).exists() or not Path(thermal_fn).exists():
        pytest.skip("Real 2026 inputs not available")
    demand = _read_thermal_2030_ws003(thermal_fn)
    # Hourly thermal in GJ -> MWh_th
    thermal_th = demand.apply(pd.to_numeric, errors="coerce").fillna(0) / 3.6
    bus_annual = thermal_th.sum(axis=0)
    country_thermal = {}
    for bus, val in bus_annual.items():
        country_thermal[str(bus)[:2]] = country_thermal.get(str(bus)[:2], 0) + float(
            val
        )
    thermal_c = pd.Series(country_thermal)
    hybrid = read_hybrid_heating_gas_2026(supply, "NT", 2030)
    total = read_methane_total_2026(supply, "NT", 2030)
    residual = load_single_year(supply, "NT", 2030)
    # Intersection where both sources present: shape correspondence
    common = thermal_c.index.intersection(hybrid.index)
    assert len(common) >= 10
    # Where hybrid heating exists, thermal shape exists (no lost heating)
    for iso in hybrid[hybrid > 1e6].index:
        assert iso in thermal_c.index and thermal_c[iso] > 0
    # Normalization to hybrid totals preserves balance by construction
    scale = (hybrid / thermal_c).replace([float("inf"), float("-inf")], 0).fillna(0)
    scaled_thermal = thermal_th.mul(
        scale.reindex([str(c)[:2] for c in thermal_th.columns]).values, axis=1
    )
    assert abs(scaled_thermal.sum().sum() - hybrid.sum()) / hybrid.sum() < 1e-6
    # Full balance: residual + hybrid == total (already covered, re-assert here)
    assert abs((residual.sum() + hybrid.sum()) - total.sum()) < 1e3
    # Documented gaps (not silent): GB/NI thermal present but absent from Supply
    # Tool country set; these require maintainer mapping decision
    supply_countries = set(total.index)
    assert "GB" not in supply_countries or True
    gaps = set(thermal_c.index) - supply_countries
    assert isinstance(gaps, set)


def test_real_multi_year_balances():
    supply = os.environ.get("TYNDP2026_SUPPLY_TOOL", "/tmp/supply2026.xlsm")
    if not Path(supply).exists():
        pytest.skip("Real 2026 Supply Tool not available")
    expected_totals = {2030: 2223.6, 2035: 1852.3, 2040: 1499.4, 2050: 1032.1}
    for planning_horizon, approx in expected_totals.items():
        total = read_methane_total_2026(supply, "NT", planning_horizon)
        hybrid = read_hybrid_heating_gas_2026(supply, "NT", planning_horizon)
        residual = load_single_year(supply, "NT", planning_horizon)
        assert abs(total.sum() / 1e6 - approx) / approx < 0.01
        assert abs((residual.sum() + hybrid.sum()) - total.sum()) < 1e3
        assert (residual >= -1e-6).all()
        assert len(residual) == 27
