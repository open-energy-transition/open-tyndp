# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""Regression tests for CBA benchmark plotting with custom projects."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pytest

matplotlib.use("Agg")


def _import_plot_module():
    """Import plot module with mocked heavy deps if needed."""
    import sys

    # Mock heavy deps only if not already available (e.g., in pixi env they are available)
    for mod in ["atlite", "pypsa", "geopandas", "xarray", "rioxarray"]:
        if mod not in sys.modules:
            sys.modules[mod] = MagicMock()
    # Mock helpers if needed
    if "scripts._helpers" not in sys.modules:
        mock_helpers = MagicMock()
        mock_helpers.configure_logging = MagicMock()
        mock_helpers.set_scenario_config = MagicMock()
        sys.modules["scripts._helpers"] = mock_helpers
    if "scripts.sb._helpers" not in sys.modules:
        mock_sb = MagicMock()
        mock_sb.add_metadata = MagicMock()
        sys.modules["scripts.sb._helpers"] = mock_sb

    from scripts.cba.plot_benchmark_indicators import (
        benchmark_range,
        plot_project_benchmarks,
        plot_summary_projects_benchmark,
    )

    return benchmark_range, plot_project_benchmarks, plot_summary_projects_benchmark


def _make_df(rows):
    return pd.DataFrame(rows)


def test_benchmark_range_missing_returns_none():
    benchmark_range, _, _ = _import_plot_module()
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
    assert result is None

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
    result2 = benchmark_range(empty, "NON_EXISTENT", source="TYNDP 2024")
    assert result2 is None


def test_benchmark_range_returns_values_when_present():
    benchmark_range, _, _ = _import_plot_module()
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


def _get_legend_texts(fig):
    texts = []
    for leg in fig.findobj(matplotlib.legend.Legend):
        texts.extend([t.get_text() for t in leg.get_texts()])
    return texts


def test_plot_project_benchmarks_custom_without_benchmark(tmp_path: Path):
    _, plot_project_benchmarks, _ = _import_plot_module()
    rows = []
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
    original_savefig = plt.Figure.savefig
    captured = {}

    def fake_savefig(self, *args, **kwargs):
        captured["fig"] = self
        return original_savefig(self, *args, **kwargs)

    plt.Figure.savefig = fake_savefig
    try:
        plot_project_benchmarks(df, out, project_label="t1500_2030", area_subtitle=None)
    finally:
        plt.Figure.savefig = original_savefig

    assert out.exists() and out.stat().st_size > 0
    fig = captured.get("fig")
    assert fig is not None
    subplot_axes = [ax for ax in fig.axes if ax.get_ylabel()]
    assert len(subplot_axes) == 3
    legends = _get_legend_texts(fig)
    assert "Open-TYNDP (mean ± min/max)" in legends
    assert "2024 TYNDP (mean ± min/max)" not in legends
    assert "Open-TYNDP low" in legends or "Open-TYNDP central" in legends
    assert "2024 TYNDP low" not in legends


def test_plot_project_benchmarks_preserves_comparison(tmp_path: Path):
    _, plot_project_benchmarks, _ = _import_plot_module()
    rows = []
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
    original_savefig = plt.Figure.savefig
    captured = {}

    def fake_savefig(self, *args, **kwargs):
        captured["fig"] = self
        return original_savefig(self, *args, **kwargs)

    plt.Figure.savefig = fake_savefig
    try:
        plot_project_benchmarks(df, out, project_label="t335_2030", area_subtitle="whole area")
    finally:
        plt.Figure.savefig = original_savefig

    assert out.exists() and out.stat().st_size > 0
    legends = _get_legend_texts(captured["fig"])
    assert "2024 TYNDP (mean ± min/max)" in legends
    assert "Open-TYNDP (mean ± min/max)" in legends
    assert "2024 TYNDP low" in legends
    assert "Open-TYNDP low" in legends


def test_plot_project_benchmarks_mixed_indicators(tmp_path: Path):
    _, plot_project_benchmarks, _ = _import_plot_module()
    rows = []
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
    for src in ["Open-TYNDP", "TYNDP 2024"]:
        for lvl, v in [("low", -2), ("central", -5), ("high", -8)]:
            rows.append(
                {
                    "source": src,
                    "indicator": "B2a_societal_cost_variation",
                    "subindex": lvl,
                    "value": v,
                    "units": "MEUR",
                    "project_id": 1501,
                    "project_code": "t1501",
                }
            )
    for src, val in [("Open-TYNDP", 5), ("TYNDP 2024", 6)]:
        rows.append(
            {
                "source": src,
                "indicator": "B3_test",
                "subindex": "mean",
                "value": val,
                "units": "GWh",
                "project_id": 1501,
                "project_code": "t1501",
            }
        )
        rows.append(
            {
                "source": src,
                "indicator": "B3_test",
                "subindex": "min",
                "value": val - 1,
                "units": "GWh",
                "project_id": 1501,
                "project_code": "t1501",
            }
        )
        rows.append(
            {
                "source": src,
                "indicator": "B3_test",
                "subindex": "max",
                "value": val + 1,
                "units": "GWh",
                "project_id": 1501,
                "project_code": "t1501",
            }
        )

    df = _make_df(rows)
    out = tmp_path / "mixed.png"
    original_savefig = plt.Figure.savefig
    captured = {}

    def fake_savefig(self, *args, **kwargs):
        captured["fig"] = self
        return original_savefig(self, *args, **kwargs)

    plt.Figure.savefig = fake_savefig
    try:
        plot_project_benchmarks(df, out)
    finally:
        plt.Figure.savefig = original_savefig

    assert out.exists()
    fig = captured["fig"]
    subplot_axes = [ax for ax in fig.axes if ax.get_ylabel()]
    assert len(subplot_axes) == 3
    legends = _get_legend_texts(fig)
    assert "Open-TYNDP (mean ± min/max)" in legends
    assert "2024 TYNDP (mean ± min/max)" in legends
    assert "Open-TYNDP central" in legends
    assert "2024 TYNDP central" in legends


def test_summary_plot_requires_both_sources(tmp_path: Path):
    _, _, plot_summary_projects_benchmark = _import_plot_module()
    rows = []
    for src in ["Open-TYNDP", "TYNDP 2024"]:
        for sub, v in [("mean", -5), ("min", -10), ("max", 0)]:
            rows.append(
                {
                    "source": src,
                    "indicator": "B1_total_system_cost_change",
                    "subindex": sub,
                    "value": v,
                    "units": "MEUR",
                    "project_id": 335,
                    "project_code": "t335",
                    "cyear": "weighted-average",
                }
            )
    for sub, v in [("mean", -7), ("min", -12), ("max", -2)]:
        rows.append(
            {
                "source": "Open-TYNDP",
                "indicator": "B1_total_system_cost_change",
                "subindex": sub,
                "value": v,
                "units": "MEUR",
                "project_id": 1500,
                "project_code": "t1500",
                "cyear": "weighted-average",
            }
        )
    df = _make_df(rows)
    out = tmp_path / "summary.png"
    original_savefig = plt.Figure.savefig
    captured = {}

    def fake_savefig(self, *args, **kwargs):
        captured["fig"] = self
        return original_savefig(self, *args, **kwargs)

    plt.Figure.savefig = fake_savefig
    try:
        plot_summary_projects_benchmark(df, out, planning_horizon="2030")
    finally:
        plt.Figure.savefig = original_savefig

    assert out.exists()
    fig = captured["fig"]
    texts = [t.get_text() for ax in fig.axes for t in ax.texts]
    assert any("n = 1" in txt for txt in texts)
    assert not any("n = 2" in txt for txt in texts)
