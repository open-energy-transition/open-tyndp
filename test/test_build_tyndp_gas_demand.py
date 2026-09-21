# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""Tests for TYNDP 2026 gas demand with Supply Tool All data and hybrid heating."""

import tempfile
from pathlib import Path

import pandas as pd
import pytest


# Import helpers after mocking dependencies
import sys, types

# Mock heavy deps if needed
for mod in ["atlite", "pypsa", "geopandas"]:
    if mod not in sys.modules:
        sys.modules[mod] = types.ModuleType(mod)

# Ensure worktree path
sys.path.insert(0, "/tmp/open-tyndp-966")

from scripts.sb.build_tyndp_gas_demand import (
    AVAILABLE_YEARS_TYNDP2026,
    _unit_factor,
    compute_thermal_ch4_annual,
    read_all_data_2026,
)


def test_available_years_include_2035_2050():
    assert 2030 in AVAILABLE_YEARS_TYNDP2026
    assert 2035 in AVAILABLE_YEARS_TYNDP2026
    assert 2040 in AVAILABLE_YEARS_TYNDP2026
    assert 2050 in AVAILABLE_YEARS_TYNDP2026
    assert len(AVAILABLE_YEARS_TYNDP2026) == 4


def test_unit_factor_twh():
    assert _unit_factor("TWh/year") == 1e6
    assert _unit_factor("TWh") == 1e6
    assert _unit_factor("GWh") == 1e3
    assert _unit_factor("MWh") == 1


def test_read_all_data_2026(tmp_path: Path):
    # Create minimal All data sheet
    data = {
        "category": ["Final demand", "Final demand", "Final demand"],
        "parameter": ["Natural gas", "Biomethane", "E-Methane"],
        "unit": ["TWh/year", "TWh/year", "TWh/year"],
        "year": [2030, 2030, 2030],
        "AT": [1.0, 0.5, 0.2],
        "DE": [10.0, 2.0, 1.0],
        "EU27": [100, 20, 10],
        "EU27.1": [100, 20, 10],
    }
    df = pd.DataFrame(data)
    fn = tmp_path / "supply.xlsx"
    df.to_excel(fn, sheet_name="All data", index=False, engine="openpyxl")

    series = read_all_data_2026(str(fn), "NT", 2030)
    assert not series.empty
    # AT: (1+0.5+0.2)=1.7 TWh -> 1.7e6 MWh
    assert abs(series["AT"] - 1.7e6) < 1
    assert abs(series["DE"] - 13e6) < 1
    # DE and AT only, EU27 ignored
    assert "EU27" not in series.index


def test_read_all_data_2026_year_filter(tmp_path: Path):
    data = {
        "category": ["Final demand", "Final demand"],
        "parameter": ["Natural gas", "Natural gas"],
        "unit": ["TWh/year", "TWh/year"],
        "year": [2030, 2040],
        "AT": [1.0, 2.0],
        "DE": [10.0, 20.0],
    }
    df = pd.DataFrame(data)
    fn = tmp_path / "supply.xlsx"
    df.to_excel(fn, sheet_name="All data", index=False, engine="openpyxl")
    s2030 = read_all_data_2026(str(fn), "NT", 2030)
    s2040 = read_all_data_2026(str(fn), "NT", 2040)
    assert s2030["AT"] == 1e6
    assert s2040["AT"] == 2e6


def test_compute_thermal_ch4_annual(tmp_path: Path):
    # Hourly thermal_ch4: 2 buses, 8760 hours
    import numpy as np

    idx = pd.date_range("2030-01-01", periods=24, freq="h")
    df = pd.DataFrame(
        {"AT00": [10.0] * 24, "DE00": [20.0] * 24, "FR00": [0] * 24}, index=idx
    )
    fn = tmp_path / "thermal_ch4.csv"
    df.to_csv(fn)

    series = compute_thermal_ch4_annual(str(fn))
    # AT: 10*24=240, DE:480, FR:0
    assert series["AT"] == 240
    assert series["DE"] == 480
    assert "FR" not in series or series["FR"] == 0


def test_subtraction_no_double_count(tmp_path: Path):
    # Integration test: supply total minus thermal = residual, energy balance
    from scripts.sb.build_tyndp_gas_demand import load_single_year

    # Create supply file
    data = {
        "category": ["Final demand"],
        "parameter": ["Natural gas"],
        "unit": ["TWh/year"],
        "year": [2030],
        "AT": [1.0],  # 1 TWh
        "DE": [2.0],  # 2 TWh
    }
    df = pd.DataFrame(data)
    fn = tmp_path / "supply.xlsx"
    df.to_excel(fn, sheet_name="All data", index=False, engine="openpyxl")

    # Create thermal file: AT 0.2 TWh, DE 0.5 TWh
    # Thermal is MW_th hourly, sum to MWh: need 0.2 TWh = 200000 MWh -> 200000/8760 ~22.83 MW avg
    # For test, use small numbers: AT 0.2e6, DE 0.5e6
    idx = pd.date_range("2030-01-01", periods=10, freq="h")
    # Use values that sum to known totals
    thermal_df = pd.DataFrame({"AT00": [20000] * 10, "DE00": [50000] * 10}, index=idx)
    thermal_fn = tmp_path / "thermal.csv"
    thermal_df.to_csv(thermal_fn)

    # Supply total: AT 1e6, DE 2e6
    # Thermal total: AT 200000, DE 500000
    # Residual should be 800000, 1500000
    result = load_single_year(str(fn), "NT", 2030, thermal_ch4_path=str(thermal_fn))
    assert abs(result["AT"] - 800000) < 1
    assert abs(result["DE"] - 1500000) < 1
    # Energy balance: residual + thermal == supply
    assert abs((result["AT"] + 200000) - 1e6) < 1
    assert abs((result["DE"] + 500000) - 2e6) < 1


def test_clipping_negative_residual(tmp_path: Path):
    from scripts.sb.build_tyndp_gas_demand import load_single_year

    data = {
        "category": ["Final demand"],
        "parameter": ["Natural gas"],
        "unit": ["TWh/year"],
        "year": [2030],
        "AT": [0.1],  # 0.1 TWh = 100000 MWh
    }
    df = pd.DataFrame(data)
    fn = tmp_path / "supply.xlsx"
    df.to_excel(fn, sheet_name="All data", index=False, engine="openpyxl")

    idx = pd.date_range("2030-01-01", periods=10, freq="h")
    thermal_df = pd.DataFrame({"AT00": [20000] * 10}, index=idx)  # 200000 MWh > supply
    thermal_fn = tmp_path / "thermal.csv"
    thermal_df.to_csv(thermal_fn)

    result = load_single_year(str(fn), "NT", 2030, thermal_ch4_path=str(thermal_fn))
    # Should clip to 0, not negative
    assert result["AT"] == 0


def test_backward_compat_fallback_2024(tmp_path: Path):
    # If All data sheet missing, should fallback to None and then 2024 logic (empty without 2024 file)
    # Here we create a file without All data sheet
    df = pd.DataFrame({"A": [1]})
    fn = tmp_path / "supply2.xlsx"
    df.to_excel(fn, sheet_name="NT+ data", index=False, engine="openpyxl")
    result = read_all_data_2026(str(fn), "NT", 2030)
    assert result is None
