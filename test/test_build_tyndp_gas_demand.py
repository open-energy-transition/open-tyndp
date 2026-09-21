# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""Tests for TYNDP 2026 gas demand with Supply Tool All data and hybrid heating."""

from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest


def _import_gas_module():
    import sys

    for mod in ["atlite", "pypsa", "geopandas", "xarray", "rioxarray", "fiona", "rasterio", "shapely"]:
        if mod not in sys.modules:
            sys.modules[mod] = MagicMock()
    if "scripts._helpers" not in sys.modules:
        mock_helpers = MagicMock()
        mock_helpers.configure_logging = MagicMock()
        mock_helpers.set_scenario_config = MagicMock()
        mock_helpers.get_snapshots = MagicMock(return_value=None)
        sys.modules["scripts._helpers"] = mock_helpers
    if "scripts.sb._helpers" not in sys.modules:
        mock_sb = MagicMock()
        # Provide dummy interpolate that does simple linear
        def dummy_interpolate(available_years, pyear, load_single_year_func, **kwargs):
            import bisect
            idx = bisect.bisect_right(sorted(available_years), pyear)
            if idx == 0:
                return load_single_year_func(pyear=available_years[0], **kwargs)
            if idx == len(available_years):
                return load_single_year_func(pyear=available_years[-1], **kwargs)
            lo = available_years[idx - 1]
            hi = available_years[idx]
            df_lo = load_single_year_func(pyear=lo, **kwargs)
            df_hi = load_single_year_func(pyear=hi, **kwargs)
            w = (pyear - lo) / (hi - lo)
            df_lo, df_hi = df_lo.align(df_hi, fill_value=0, join="outer")
            return df_lo * (1 - w) + df_hi * w

        mock_sb.interpolate_demand = dummy_interpolate
        mock_sb.safe_pyear = MagicMock(return_value=2030)
        sys.modules["scripts.sb._helpers"] = mock_sb

    from scripts.sb.build_tyndp_gas_demand import (
        AVAILABLE_YEARS_TYNDP2026,
        DEFAULT_GAS_BOILER_EFFICIENCY,
        _unit_factor,
        compute_thermal_ch4_annual,
        read_hybrid_heating_gas_2026,
        read_methane_total_2026,
    )

    return (
        AVAILABLE_YEARS_TYNDP2026,
        DEFAULT_GAS_BOILER_EFFICIENCY,
        _unit_factor,
        compute_thermal_ch4_annual,
        read_hybrid_heating_gas_2026,
        read_methane_total_2026,
    )


def test_available_years_include_2035_2050():
    AVAILABLE_YEARS_TYNDP2026, _, _, _, _, _ = _import_gas_module()
    assert 2030 in AVAILABLE_YEARS_TYNDP2026
    assert 2035 in AVAILABLE_YEARS_TYNDP2026
    assert 2040 in AVAILABLE_YEARS_TYNDP2026
    assert 2050 in AVAILABLE_YEARS_TYNDP2026


def test_unit_factor():
    _, _, _unit_factor, _, _, _ = _import_gas_module()
    assert _unit_factor("TWh/year") == 1e6
    assert _unit_factor("TWh") == 1e6
    assert _unit_factor("GWh") == 1e3
    assert _unit_factor("MWh") == 1


def test_read_methane_total_2026(tmp_path):
    _, _, _, _, _, read_methane_total_2026 = _import_gas_module()
    df = pd.DataFrame(
        {
            "ETM": ["Methane"],
            "Parameter": ["Total Energy Demand"],
            "Unit": ["TWh/year"],
            "Year": [2030],
            "AT": [1.7],
            "DE": [13.0],
            "EU27": [100.0],
            "EU27.1": [100.0],
        }
    )
    fn = tmp_path / "supply.xlsx"
    df.to_excel(fn, sheet_name="All data", index=False, engine="openpyxl")
    series = read_methane_total_2026(str(fn), "NT", 2030)
    assert not series.empty
    assert abs(series["AT"] - 1.7e6) < 1
    assert abs(series["DE"] - 13e6) < 1
    assert "EU27" not in series.index


def test_read_hybrid_heating_2026(tmp_path):
    _, _, _, _, read_hybrid_heating_gas_2026, _ = _import_gas_module()
    df = pd.DataFrame(
        {
            "ETM": ["Methane"],
            "Parameter": ["Methane for hybrid heating"],
            "Unit": ["TWh/year"],
            "Year": [2030],
            "AT": [0.5],
            "DE": [1.35],
            "EU27": [33.3],
        }
    )
    fn = tmp_path / "supply.xlsx"
    df.to_excel(fn, sheet_name="All data", index=False, engine="openpyxl")
    series = read_hybrid_heating_gas_2026(str(fn), "NT", 2030)
    assert abs(series["AT"] - 0.5e6) < 1
    assert abs(series["DE"] - 1.35e6) < 1


def test_read_all_data_year_filter(tmp_path):
    _, _, _, _, _, read_methane_total_2026 = _import_gas_module()
    df = pd.DataFrame(
        {
            "ETM": ["Methane", "Methane"],
            "Parameter": ["Total Energy Demand", "Total Energy Demand"],
            "Unit": ["TWh/year", "TWh/year"],
            "Year": [2030, 2040],
            "AT": [1.0, 2.0],
            "DE": [10.0, 20.0],
        }
    )
    fn = tmp_path / "supply.xlsx"
    df.to_excel(fn, sheet_name="All data", index=False, engine="openpyxl")
    s2030 = read_methane_total_2026(str(fn), "NT", 2030)
    s2040 = read_methane_total_2026(str(fn), "NT", 2040)
    assert s2030["AT"] == 1e6
    assert s2040["AT"] == 2e6


def test_compute_thermal_ch4_annual(tmp_path):
    _, _, _, compute_thermal_ch4_annual, _, _ = _import_gas_module()
    idx = pd.date_range("2030-01-01", periods=24, freq="h")
    df = pd.DataFrame({"AT00": [10.0] * 24, "DE00": [20.0] * 24}, index=idx)
    fn = tmp_path / "thermal.csv"
    df.to_csv(fn)
    series = compute_thermal_ch4_annual(str(fn), efficiency=1.0)
    assert series["AT"] == 240
    assert series["DE"] == 480
    series_eff = compute_thermal_ch4_annual(str(fn), efficiency=0.9)
    assert abs(series_eff["AT"] - 240 / 0.9) < 1e-6


def test_subtraction_uses_supply_tool_hybrid_not_thermal(tmp_path):
    from scripts.sb.build_tyndp_gas_demand import load_single_year

    df = pd.DataFrame(
        {
            "ETM": ["Methane", "Methane"],
            "Parameter": ["Total Energy Demand", "Methane for hybrid heating"],
            "Unit": ["TWh/year", "TWh/year"],
            "Year": [2030, 2030],
            "AT": [1.0, 0.2],
            "DE": [2.0, 0.5],
        }
    )
    fn = tmp_path / "supply.xlsx"
    df.to_excel(fn, sheet_name="All data", index=False, engine="openpyxl")
    idx = pd.date_range("2030-01-01", periods=10, freq="h")
    thermal_df = pd.DataFrame({"AT00": [1000] * 10, "DE00": [1000] * 10}, index=idx)
    thermal_fn = tmp_path / "thermal.csv"
    thermal_df.to_csv(thermal_fn)
    result = load_single_year(str(fn), "NT", 2030, thermal_ch4_path=str(thermal_fn))
    assert abs(result["AT"] - 800000) < 1
    assert abs(result["DE"] - 1500000) < 1


def test_negative_residual_raises(tmp_path):
    from scripts.sb.build_tyndp_gas_demand import load_single_year

    df = pd.DataFrame(
        {
            "ETM": ["Methane", "Methane"],
            "Parameter": ["Total Energy Demand", "Methane for hybrid heating"],
            "Unit": ["TWh/year", "TWh/year"],
            "Year": [2030, 2030],
            "AT": [0.1, 0.2],
        }
    )
    fn = tmp_path / "supply.xlsx"
    df.to_excel(fn, sheet_name="All data", index=False, engine="openpyxl")
    idx = pd.date_range("2030-01-01", periods=10, freq="h")
    thermal_df = pd.DataFrame({"AT00": [100] * 10}, index=idx)
    thermal_fn = tmp_path / "thermal.csv"
    thermal_df.to_csv(thermal_fn)
    with pytest.raises(ValueError, match="Negative residual"):
        load_single_year(str(fn), "NT", 2030, thermal_ch4_path=str(thermal_fn))


def test_missing_thermal_for_2026_still_subtracts_hybrid(tmp_path):
    from scripts.sb.build_tyndp_gas_demand import load_single_year

    df = pd.DataFrame(
        {
            "ETM": ["Methane", "Methane"],
            "Parameter": ["Total Energy Demand", "Methane for hybrid heating"],
            "Unit": ["TWh/year", "TWh/year"],
            "Year": [2030, 2030],
            "AT": [1.0, 0.2],
        }
    )
    fn = tmp_path / "supply.xlsx"
    df.to_excel(fn, sheet_name="All data", index=False, engine="openpyxl")
    result = load_single_year(str(fn), "NT", 2030, thermal_ch4_path=None)
    assert abs(result["AT"] - 800000) < 1


def test_snakemake_rule_has_thermal_input():
    rule_file = Path(__file__).parent.parent / "rules/sb.smk"
    if not rule_file.exists():
        pytest.skip("sb.smk not found")
    content = rule_file.read_text()
    assert "thermal_ch4" in content
    assert "build_tyndp_demand" in content
    assert "demand_tyndp_thermal_ch4" in content


def test_real_supply_tool_file():
    real = Path("/tmp/supply2026.xlsm")
    if not real.exists():
        pytest.skip("Real Supply Tool not downloaded")
    _, _, _, _, read_hybrid_heating_gas_2026, read_methane_total_2026 = _import_gas_module()
    total = read_methane_total_2026(str(real), "NT", 2030)
    hybrid = read_hybrid_heating_gas_2026(str(real), "NT", 2030)
    assert not total.empty
    assert total["DE"] > 400e6
    assert hybrid["DE"] > 1e6
    residual = total["DE"] - hybrid["DE"]
    assert residual > 400e6


def test_year_specific_interpolation(tmp_path):
    from scripts.sb.build_tyndp_gas_demand import load_gas_demand

    # Create supply with 2030 and 2035, request 2032 (interpolated)
    # 2030 residual 0.8, 2035 residual 1.0, 2032 weight 0.4 => 0.88
    df = pd.DataFrame(
        {
            "ETM": ["Methane", "Methane", "Methane", "Methane"],
            "Parameter": ["Total Energy Demand", "Methane for hybrid heating", "Total Energy Demand", "Methane for hybrid heating"],
            "Unit": ["TWh/year", "TWh/year", "TWh/year", "TWh/year"],
            "Year": [2030, 2030, 2035, 2035],
            "AT": [1.0, 0.2, 1.5, 0.5],
        }
    )
    fn = tmp_path / "supply.xlsx"
    df.to_excel(fn, sheet_name="All data", index=False, engine="openpyxl")

    # Request 2032, between 2030 and 2035 (weight 0.4)
    # 2030 residual 0.8, 2035 residual 1.0 => 0.8*0.6 + 1.0*0.4 = 0.88
    result = load_gas_demand(str(fn), "NT", 2032, thermal_ch4_path=None)
    assert abs(result["AT"] - 880000) < 1
