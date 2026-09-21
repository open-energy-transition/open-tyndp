# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""Regression tests for CBA benchmark plotting with custom projects."""

import tempfile
from pathlib import Path

import pandas as pd
import pytest

from scripts.cba.plot_benchmark_indicators import (
    benchmark_range,
    plot_project_benchmarks,
)


def _make_df(rows):
    """Helper to build indicator DataFrame."""
    return pd.DataFrame(rows)


def test_benchmark_range_missing_returns_tuple():
    df = _make_df(
        [
            {
                "source": "Open-TYNDP",
                "indicator": "B1_total_system_cost_change",
                "subindex": "mean",
                "value": 1.0,
                "project_id": 999,
            }
        ]
    )
    result = benchmark_range(df, "B1_total_system_cost_change", source="TYNDP 2024")
    assert result == (None, None, None), f"expected (None,None,None) got {result}"

    empty = _make_df(
        [
            {
                "source": "TYNDP 2024",
                "indicator": "B1_total_system_cost_change",
                "subindex": "mean",
                "value": 1.0,
                "project_id": 999,
            }
        ]
    )
    # Query for non-existent indicator should also return tuple
    result2 = benchmark_range(empty, "NON_EXISTENT", source="TYNDP 2024")
    assert result2 == (None, None, None)


def test_benchmark_range_returns_values_when_present():
    df = _make_df(
        [
            {"source": "TYNDP 2024", "indicator": "B1_total_system_cost_change", "subindex": "min", "value": -10, "project_id": 335},
            {"source": "TYNDP 2024", "indicator": "B1_total_system_cost_change", "subindex": "mean", "value": -5, "project_id": 335},
            {"source": "TYNDP 2024", "indicator": "B1_total_system_cost_change", "subindex": "max", "value": 0, "project_id": 335},
        ]
    )
    mn, mean, mx = benchmark_range(df, "B1_total_system_cost_change", source="TYNDP 2024")
    assert mn == -10
    assert mean == -5
    assert mx == 0


def test_plot_project_benchmarks_custom_without_benchmark(tmp_path: Path):
    """Custom project without TYNDP benchmark should plot Open-TYNDP markers."""
    rows = []
    # Custom project 1500 with all indicators but only Open-TYNDP source
    indicators = [
        ("B1_total_system_cost_change", "mean", -123.4, "MEUR"),
        ("B1_total_system_cost_change", "min", -150, "MEUR"),
        ("B1_total_system_cost_change", "max", -100, "MEUR"),
        ("B2a_societal_cost_variation", "low", -10, "MEUR"),
        ("B2a_societal_cost_variation", "central", -5, "MEUR"),
        ("B2a_societal_cost_variation", "high", -1, "MEUR"),
        ("B3_res Curtailment", "mean", 5, "GWh"),
        ("B3_res Curtailment", "min", 2, "GWh"),
        ("B3_res Curtailment", "max", 8, "GWh"),
    ]
    for ind, sub, val, unit in indicators:
        rows.append(
            {
                "source": "Open-TYNDP",
                "indicator": ind,
                "subindex": sub,
                "value": val,
                "units": unit,
                "project_id": 1500,
                "project_code": "t1500",
            }
        )
    df = _make_df(rows)
    out = tmp_path / "custom.png"
    plot_project_benchmarks(df, out, project_label="t1500_2030", area_subtitle=None)
    assert out.exists(), "Plot not created for custom project without benchmark"
    assert out.stat().st_size > 0


def test_plot_project_benchmarks_preserves_comparison(tmp_path: Path):
    """Official project with both sources should plot both markers."""
    rows = []
    # B1 with both sources
    for src, val in [("Open-TYNDP", -5), ("TYNDP 2024", -6)]:
        for sub, v in [("mean", val), ("min", val - 5), ("max", val + 5)]:
            rows.append(
                {
                    "source": src,
                    "indicator": "B1_total_system_cost_change",
                    "subindex": sub,
                    "value": v,
                    "units": "MEUR",
                    "project_id": 335,
                    "project_code": "t335",
                }
            )
    # B2a with both
    for src in ["Open-TYNDP", "TYNDP 2024"]:
        for lvl, v in [("low", -2), ("central", -5), ("high", -8)]:
            rows.append(
                {
                    "source": src,
                    "indicator": "B2a_societal_cost_variation",
                    "subindex": lvl,
                    "value": v,
                    "units": "MEUR",
                    "project_id": 335,
                    "project_code": "t335",
                }
            )
    df = _make_df(rows)
    out = tmp_path / "official.png"
    plot_project_benchmarks(df, out, project_label="t335_2030", area_subtitle="whole area")
    assert out.exists()
    assert out.stat().st_size > 0


def test_plot_project_benchmarks_mixed_indicators(tmp_path: Path):
    """Mixed: B1 without benchmark, B2a with benchmark should still plot both."""
    rows = []
    # B1 only Open-TYNDP
    for sub, v in [("mean", -10), ("min", -15), ("max", -5)]:
        rows.append(
            {
                "source": "Open-TYNDP",
                "indicator": "B1_total_system_cost_change",
                "subindex": sub,
                "value": v,
                "units": "MEUR",
                "project_id": 1501,
                "project_code": "t1501",
            }
        )
    # B2a both sources
    for src in ["Open-TYNDP", "TYNDP 2024"]:
        rows.append(
            {
                "source": src,
                "indicator": "B2a_societal_cost_variation",
                "subindex": "central",
                "value": -5,
                "units": "MEUR",
                "project_id": 1501,
                "project_code": "t1501",
            }
        )
    df = _make_df(rows)
    out = tmp_path / "mixed.png"
    plot_project_benchmarks(df, out)
    assert out.exists()
