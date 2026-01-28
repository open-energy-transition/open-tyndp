# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Create benchmark plots for CBA indicators.

This script reads a indicators CSV file that
includes model rows and TYNDP benchmark rows,
and generates plots comparing model values
to the 2024 TYNDP values.
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

def select_model_value(df: pd.DataFrame, indicator: str) -> float | None:
    """Pick a representative model value for a given indicator."""
    model = df[(df["source"] == "model") & (df["indicator"] == indicator)]
    if model.empty:
        return None

    for subindex in ["mean", "central", ""]:
        subset = model[model["subindex"] == subindex]
        if not subset.empty:
            return subset["value"].mean()

    return model["value"].mean()


def select_value_by_subindex(
    df: pd.DataFrame, indicator: str, source: str, subindex: str
) -> float | None:
    """Return the mean value for a given indicator/subindex/source."""
    subset = df[
        (df["source"] == source)
        & (df["indicator"] == indicator)
        & (df["subindex"] == subindex)
    ]
    if subset.empty:
        return None
    return float(subset["value"].mean())


def benchmark_range(
    df: pd.DataFrame, indicator: str
) -> tuple[float, float, float] | None:
    """Return (min, mean, max) range for a benchmark indicator."""
    benchmark = df[
        (df["source"] == "2024 tyndp") & (df["indicator"] == indicator)
    ].copy()
    if benchmark.empty:
        return None

    if indicator == "B2a_societal_cost_variation":
        benchmark = benchmark[benchmark["subindex"].isin(["low", "central", "high"])]
    else:
        benchmark = benchmark[
            benchmark["subindex"].isin(["min", "mean", "max", "explicit"])
        ]

    benchmark = benchmark.assign(
        subindex=benchmark["subindex"].replace({"explicit": "mean"})
    )
    pivot = benchmark.pivot_table(
        index="project_id", columns="subindex", values="value", aggfunc="mean"
    )
    if pivot.empty:
        return None

    mean = pivot.get("mean")
    if mean is None:
        if "min" in pivot.columns:
            mean = pivot["min"]
        elif "max" in pivot.columns:
            mean = pivot["max"]
    if mean is None or mean.empty:
        return None

    mean_val = float(mean.iloc[0])
    min_val = float(pivot.get("min", mean).iloc[0])
    max_val = float(pivot.get("max", mean).iloc[0])
    return min_val, mean_val, max_val


def format_project_label(project_id: int, planning_horizon: str | None) -> str:
    """Format plot title label from project id and planning horizon."""
    if planning_horizon:
        return f"t{project_id}_{planning_horizon}"
    return f"t{project_id}"


def plot_project_benchmarks(
    df: pd.DataFrame, output_path: Path, project_label: str | None = None
) -> None:
    """Plot one subplot per indicator with its own y-axis and legend."""
    indicators = sorted(df["indicator"].dropna().unique())
    if not indicators:
        logger.info("No benchmark indicators available to plot")
        return

    plot_items = []
    for indicator in indicators:
        if indicator == "B2a_societal_cost_variation":
            levels = ["low", "central", "high"]
            has_any = any(
                select_value_by_subindex(df, indicator, src, lvl) is not None
                for src in ["model", "2024 tyndp"]
                for lvl in levels
            )
            if has_any:
                plot_items.append(indicator)
            continue

        model_val = select_model_value(df, indicator)
        bench = benchmark_range(df, indicator)
        if model_val is None or bench is None:
            continue
        plot_items.append(indicator)

    if not plot_items:
        logger.info("No complete benchmark data to plot")
        return

    fig, axes = plt.subplots(
        nrows=1,
        ncols=len(plot_items),
        figsize=(max(1.7 * len(plot_items), 10), 4.2),
        squeeze=False,
        sharey=False,
    )

    legend_handles = []
    legend_labels = []

    for ax, indicator in zip(axes[0], plot_items):
        ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.6)

        if indicator == "B2a_societal_cost_variation":
            levels = ["low", "central", "high"]
            level_colors = {
                "low": "tab:orange",
                "central": "tab:green",
                "high": "tab:red",
            }
            offsets = {"low": -0.1, "central": 0.0, "high": 0.1}
            x = 0.0
            for level in levels:
                model_val = select_value_by_subindex(df, indicator, "model", level)
                bench_val = select_value_by_subindex(df, indicator, "2024 tyndp", level)
                if bench_val is not None:
                    ax.scatter(
                        [x + offsets[level]],
                        [bench_val],
                        color=level_colors[level],
                        marker="x",
                        zorder=3,
                    )
                    label = f"2024 TYNDP {level}"
                    if label not in legend_labels:
                        legend_handles.append(
                            Line2D(
                                [0],
                                [0],
                                marker="x",
                                color=level_colors[level],
                                linestyle="None",
                            )
                        )
                        legend_labels.append(label)
                if model_val is not None:
                    ax.scatter(
                        [x + offsets[level]],
                        [model_val],
                        color=level_colors[level],
                        marker="o",
                        zorder=4,
                    )
                    label = f"Model {level}"
                    if label not in legend_labels:
                        legend_handles.append(
                            Line2D(
                                [0],
                                [0],
                                marker="o",
                                color=level_colors[level],
                                linestyle="None",
                            )
                        )
                        legend_labels.append(label)
            ax.set_xticks([])
        else:
            model_val = select_model_value(df, indicator)
            min_val, mean_val, max_val = benchmark_range(df, indicator)
            lower = abs(mean_val - min_val)
            upper = abs(max_val - mean_val)

            ax.errorbar(
                [0],
                [mean_val],
                yerr=[[lower], [upper]],
                fmt="o",
                color="gray",
                ecolor="lightgray",
                capsize=3,
            )
            label = "2024 TYNDP (mean ± min/max)"
            if label not in legend_labels:
                legend_handles.append(
                    Line2D([0], [0], marker="o", color="gray", linestyle="None")
                )
                legend_labels.append(label)
            ax.scatter(
                [0],
                [model_val],
                color="tab:blue",
                zorder=3,
            )
            if "Model" not in legend_labels:
                legend_handles.append(
                    Line2D([0], [0], marker="o", color="tab:blue", linestyle="None")
                )
                legend_labels.append("Model")
            ax.set_xticks([])

        ylim = ax.get_ylim()
        values = [v for v in ylim if v != 0]
        if values:
            abs_max = max(abs(v) for v in values)
            abs_min = min(abs(v) for v in values)
            if abs_min > 0 and abs_max / abs_min >= 1e3:
                ax.set_yscale("symlog", linthresh=max(abs_min, 1.0))

        units = df.loc[
            (df["indicator"] == indicator) & (df["source"] == "model"), "units"
        ].dropna()
        unit_label = units.iloc[0] if not units.empty else ""
        title = f"{indicator} ({unit_label})" if unit_label else indicator
        ax.set_ylabel(title)

    if project_label:
        fig.suptitle(f"CBA indicator benchmark ({project_label})")
    b2a_order = [
        "2024 TYNDP low",
        "2024 TYNDP central",
        "2024 TYNDP high",
        "Model low",
        "Model central",
        "Model high",
    ]
    b2a_items = [
        (h, l) for h, l in zip(legend_handles, legend_labels) if l in b2a_order
    ]
    other_items = [
        (h, l) for h, l in zip(legend_handles, legend_labels) if l not in b2a_order
    ]

    if b2a_items:
        ordered_b2a = [(h, l) for l in b2a_order for h, l2 in b2a_items if l2 == l]
        b2a_handles = [h for h, _ in ordered_b2a]
        b2a_labels = [l for _, l in ordered_b2a]
        fig.legend(
            b2a_handles,
            b2a_labels,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.02),
            ncol=3,
            frameon=False,
        )

    if other_items:
        other_handles = [h for h, _ in other_items]
        other_labels = [l for _, l in other_items]
        fig.legend(
            other_handles,
            other_labels,
            loc="lower center",
            bbox_to_anchor=(0.12, -0.02),
            ncol=1,
            frameon=False,
        )
    fig.tight_layout(rect=[0, 0.12, 1, 0.96])
    fig.savefig(output_path, dpi=400)
    plt.close(fig)


def create_plots(indicators_file, output_path, planning_horizon=None):
    """Create benchmark plots from a per-project or collected indicators file."""
    output_path = Path(output_path)
    output_dir = output_path if output_path.suffix == "" else output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(indicators_file)
    if df.empty:
        logger.warning("No indicators data to plot")
        return

    if "source" not in df.columns:
        df["source"] = "model"

    project_ids = df.loc[df["source"] == "model", "project_id"].dropna().unique()
    if output_path.suffix:
        if len(project_ids) == 0:
            logger.info("No model projects found in indicators file")
            return
        project_id = int(project_ids[0])
        project_label = format_project_label(project_id, planning_horizon)
        project_df = df[df["project_id"] == project_id].copy()
        plot_project_benchmarks(project_df, output_path, project_label)
    else:
        for project_id in project_ids:
            project_df = df[df["project_id"] == project_id].copy()
            project_id = int(project_id)
            suffix = f"_{planning_horizon}" if planning_horizon else ""
            output_file = output_dir / f"project_{project_id}{suffix}.png"
            project_label = format_project_label(project_id, planning_horizon)
            plot_project_benchmarks(project_df, output_file, project_label)

    logger.info("Benchmark plots saved to %s", output_dir)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("plot_benchmark_indicators")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    planning_horizon = snakemake.wildcards.get("planning_horizons")
    output_target = snakemake.output.get("plot_file") or snakemake.output.plot_dir
    create_plots(snakemake.input.indicators, output_target, planning_horizon)
