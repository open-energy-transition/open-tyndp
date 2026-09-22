# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""Tests for TYNDP 2026 gas demand (All data row selection and hybrid subtraction)."""

import os
from pathlib import Path

import pandas as pd
import pytest

from scripts.sb.build_tyndp_gas_demand import (
    AVAILABLE_YEARS_TYNDP2026,
    _unit_factor,
    load_gas_demand,
    load_single_year,
    read_hybrid_heating_gas_2026,
    read_methane_total_2026,
)


def test_available_years():
    assert AVAILABLE_YEARS_TYNDP2026 == [2030, 2035, 2040, 2050]


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
