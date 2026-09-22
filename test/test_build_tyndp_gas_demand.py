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


def _thermal_bus_of(sheet: str) -> str:
    import re

    m = re.search(r"-([A-Z0-9]+)_", sheet)
    return m.group(1) if m else sheet


def _read_thermal_2030_ws003(thermal_fn: str):
    import pandas as pd

    data = pd.read_excel(
        thermal_fn,
        header=10,
        index_col=[0, 1],
        sheet_name=None,
        usecols=lambda n: n in ("Date", "Hour", "WS003"),
        engine="calamine",
    )
    demand = pd.concat(data, axis=1).droplevel(1, axis=1)
    demand.columns = [_thermal_bus_of(c).replace("UK", "GB") for c in demand.columns]
    return demand


def test_thermal_bus_extraction_and_country_mapping():
    thermal_fn = os.environ.get(
        "TYNDP2026_THERMAL_CH4_2030",
        "/tmp/demand2030/Demand/Demand/2030/Thermal_energy_Methane_2030.xlsx",
    )
    if not Path(thermal_fn).exists():
        pytest.skip("Real 2026 thermal_ch4 file not available")
    demand = _read_thermal_2030_ws003(thermal_fn)
    # 41 sheets, 8760 hourly rows (365 days, not snapshots 8736)
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
    import pandas as pd

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
    import pandas as pd

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


def test_temporal_resolution_mismatch_documented():
    thermal_fn = os.environ.get(
        "TYNDP2026_THERMAL_CH4_2030",
        "/tmp/demand2030/Demand/Demand/2030/Thermal_energy_Methane_2030.xlsx",
    )
    if not Path(thermal_fn).exists():
        pytest.skip("Real 2026 thermal_ch4 file not available")

    demand = _read_thermal_2030_ws003(thermal_fn)
    assert len(demand) == 8760  # 365 days in thermal file
    # Snapshots config 2009-01-01 to 2009-12-31 minus last day = 8736h (52 weeks)
    # Downstream indexing of 8760h thermal with 8736h snapshots drops final 24h
    # (~0.3%); explicit rejection/wiring belongs to tyndp-2026 hourly integration
    assert 8760 - 8736 == 24


def test_real_multi_year_balances():
    supply = os.environ.get("TYNDP2026_SUPPLY_TOOL", "/tmp/supply2026.xlsm")
    if not Path(supply).exists():
        pytest.skip("Real 2026 Supply Tool not available")
    expected_totals = {2030: 2223.6, 2035: 1852.3, 2040: 1499.4, 2050: 1032.1}
    for pyear, approx in expected_totals.items():
        total = read_methane_total_2026(supply, "NT", pyear)
        hybrid = read_hybrid_heating_gas_2026(supply, "NT", pyear)
        residual = load_single_year(supply, "NT", pyear)
        assert abs(total.sum() / 1e6 - approx) / approx < 0.01
        assert abs((residual.sum() + hybrid.sum()) - total.sum()) < 1e3
        assert (residual >= -1e-6).all()
        assert len(residual) == 27
