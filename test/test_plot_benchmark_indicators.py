# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""Regression tests for CBA benchmark plotting with custom projects."""

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

matplotlib.use("Agg")

from scripts.cba.plot_benchmark_indicators import (
    benchmark_range,
    plot_project_benchmarks,
    plot_summary_projects_benchmark,
)


def _make_df(rows):
    return pd.DataFrame(rows)


def _get_legend_texts(fig):
    texts = []
    for leg in fig.findobj(matplotlib.legend.Legend):
        texts.extend([t.get_text() for t in leg.get_texts()])
    return texts


def _capture_fig(monkeypatch):
    captured = {}

    def fake_savefig(self, *args, **kwargs):
        captured["fig"] = self

    monkeypatch.setattr(plt.Figure, "savefig", fake_savefig)
    return captured


def test_benchmark_range_missing_returns_none():
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
    assert (
        benchmark_range(df, "B1_total_system_cost_change", source="TYNDP 2024") is None
    )

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
    assert benchmark_range(empty, "NON_EXISTENT", source="TYNDP 2024") is None


def test_benchmark_range_returns_values_when_present():
    df = _make_df(
        [
            {
                "source": "TYNDP 2024",
                "indicator": "B1_total_system_cost_change",
                "subindex": "min",
                "value": -10,
                "project_id": 335,
            },
            {
                "source": "TYNDP 2024",
                "indicator": "B1_total_system_cost_change",
                "subindex": "mean",
                "value": -5,
                "project_id": 335,
            },
            {
                "source": "TYNDP 2024",
                "indicator": "B1_total_system_cost_change",
                "subindex": "max",
                "value": 0,
                "project_id": 335,
            },
        ]
    )
    mn, mean, mx = benchmark_range(
        df, "B1_total_system_cost_change", source="TYNDP 2024"
    )
    assert (mn, mean, mx) == (-10, -5, 0)


def test_plot_project_benchmarks_custom_without_benchmark(tmp_path: Path, monkeypatch):
    rows = [
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
    df = _make_df(
        [
            {
                "source": "Open-TYNDP",
                "indicator": ind,
                "subindex": sub,
                "value": val,
                "units": unit,
                "project_id": 1500,
                "project_code": "t1500",
            }
            for ind, sub, val, unit in rows
        ]
    )
    out = tmp_path / "custom.png"
    captured = _capture_fig(monkeypatch)
    plot_project_benchmarks(df, out, project_label="t1500_2030", area_subtitle=None)

    fig = captured["fig"]
    subplot_axes = [ax for ax in fig.axes if ax.get_ylabel()]
    assert len(subplot_axes) == 3

    # Verify actual plotted Open-TYNDP values, no TYNDP markers fabricated
    by_indicator = {ax.get_ylabel().split(" ")[0]: ax for ax in subplot_axes}
    b1_ax = by_indicator["B1_total_system_cost_change"]
    # Open-TYNDP errorbar at x=0.1 with y=-123.4 should exist
    b1_ys = [line.get_ydata()[0] for line in b1_ax.lines if len(line.get_ydata())]
    assert any(abs(y - (-123.4)) < 1e-6 for y in b1_ys)
    # No TYNDP errorbar at x=-0.1
    b1_xs = [line.get_xdata()[0] for line in b1_ax.lines if len(line.get_xdata())]
    assert not any(abs(x - (-0.1)) < 1e-6 for x in b1_xs)

    legends = _get_legend_texts(fig)
    assert "Open-TYNDP (mean ± min/max)" in legends
    assert "2024 TYNDP (mean ± min/max)" not in legends
    assert "2024 TYNDP low" not in legends


def test_plot_project_benchmarks_preserves_comparison(tmp_path: Path, monkeypatch):
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
    captured = _capture_fig(monkeypatch)
    plot_project_benchmarks(
        df, out, project_label="t335_2030", area_subtitle="whole area"
    )

    fig = captured["fig"]
    by_indicator = {
        ax.get_ylabel().split(" ")[0]: ax for ax in fig.axes if ax.get_ylabel()
    }
    b1_ax = by_indicator["B1_total_system_cost_change"]
    b1_points = [
        (line.get_xdata()[0], line.get_ydata()[0])
        for line in b1_ax.lines
        if len(line.get_xdata())
    ]
    # Both markers present with correct values
    assert any(abs(x - (-0.1)) < 1e-6 and abs(y - (-6)) < 1e-6 for x, y in b1_points)
    assert any(abs(x - 0.1) < 1e-6 and abs(y - (-5)) < 1e-6 for x, y in b1_points)

    legends = _get_legend_texts(fig)
    assert "2024 TYNDP (mean ± min/max)" in legends
    assert "Open-TYNDP (mean ± min/max)" in legends


def test_plot_project_benchmarks_mixed_indicators(tmp_path: Path, monkeypatch):
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
        for sub, delta in [("mean", 0), ("min", -1), ("max", 1)]:
            rows.append(
                {
                    "source": src,
                    "indicator": "B3_test",
                    "subindex": sub,
                    "value": val + delta,
                    "units": "GWh",
                    "project_id": 1501,
                    "project_code": "t1501",
                }
            )

    df = _make_df(rows)
    out = tmp_path / "mixed.png"
    captured = _capture_fig(monkeypatch)
    plot_project_benchmarks(df, out)

    fig = captured["fig"]
    subplot_axes = [ax for ax in fig.axes if ax.get_ylabel()]
    assert len(subplot_axes) == 3
    by_indicator = {ax.get_ylabel().split(" ")[0]: ax for ax in subplot_axes}
    # B1 without benchmark: only Open marker
    b1_xs = [
        line.get_xdata()[0]
        for line in by_indicator["B1_total_system_cost_change"].lines
        if len(line.get_xdata())
    ]
    assert not any(abs(x - (-0.1)) < 1e-6 for x in b1_xs)
    # B3 with benchmark: both markers
    b3_points = [
        (line.get_xdata()[0], line.get_ydata()[0])
        for line in by_indicator["B3_test"].lines
        if len(line.get_xdata())
    ]
    assert any(abs(x - (-0.1)) < 1e-6 for x, _ in b3_points)
    assert any(abs(x - 0.1) < 1e-6 for x, _ in b3_points)


def test_summary_plot_requires_both_sources(tmp_path: Path, monkeypatch):
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
    captured = _capture_fig(monkeypatch)
    plot_summary_projects_benchmark(df, out, planning_horizon="2030")

    fig = captured["fig"]
    texts = [t.get_text() for ax in fig.axes for t in ax.texts]
    assert any("n = 1" in txt for txt in texts)
    assert not any("n = 2" in txt for txt in texts)
