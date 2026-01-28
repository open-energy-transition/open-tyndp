# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Create benchmark plots for CBA indicators.

This script reads a indicators CSV file (single project or collected) that
includes model rows and TYNDP benchmark rows, and generates plots comparing
model values (dot) to TYNDP min/mean/max bands.
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

MODEL_SUBINDEX_PREFERENCE = {
    "B2a_societal_cost_variation": ["central", "mean", ""],
}


def select_model_value(df: pd.DataFrame, indicator: str) -> float | None:
    """Pick a representative model value for a given indicator."""
    model = df[(df["source"] == "model") & (df["indicator"] == indicator)]
    if model.empty:
        return None

    if indicator.startswith("B4"):
        candidates = ["mean", ""]
    else:
        candidates = MODEL_SUBINDEX_PREFERENCE.get(indicator, ["", "mean", "central"])

    for subindex in candidates:
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
    """Plot benchmark ranges and model values for a single project."""
    indicators = sorted(df["indicator"].dropna().unique())
    if not indicators:
        logger.info("No benchmark indicators available to plot")
        return

    labels = []
    model_values = []
    means = []
    mins = []
    maxs = []
    b2a_special = None

    for indicator in indicators:
        if indicator == "B2a_societal_cost_variation":
            b2a_special = indicator
            continue

        model_val = select_model_value(df, indicator)
        bench = benchmark_range(df, indicator)
        if model_val is None or bench is None:
            continue
        min_val, mean_val, max_val = bench
        units = df.loc[
            (df["indicator"] == indicator) & (df["source"] == "model"), "units"
        ].dropna()
        unit_label = units.iloc[0] if not units.empty else ""
        labels.append(f"{indicator} ({unit_label})" if unit_label else indicator)
        model_values.append(model_val)
        means.append(mean_val)
        mins.append(min_val)
        maxs.append(max_val)

    if not labels:
        logger.info("No complete benchmark data to plot")
        return

    x = range(len(labels))
    mean = pd.Series(means).to_numpy()
    ymin = pd.Series(mins).to_numpy()
    ymax = pd.Series(maxs).to_numpy()
    lower = abs(mean - ymin)
    upper = abs(ymax - mean)

    fig, ax = plt.subplots(figsize=(max(10, len(labels) * 0.6), 5))
    # Switch to symlog if values span multiple orders of magnitude.
    values = list(mean) + list(ymin) + list(ymax) + list(model_values)
    nonzero = [abs(v) for v in values if v not in (0, None) and not pd.isna(v)]
    if nonzero:
        abs_max = max(nonzero)
        abs_min = min(nonzero)
        if abs_min > 0 and abs_max / abs_min >= 1e3:
            ax.set_yscale("symlog", linthresh=max(abs_min, 1.0))
    ax.errorbar(
        x,
        mean,
        yerr=[lower, upper],
        fmt="o",
        color="gray",
        ecolor="lightgray",
        capsize=3,
        label="2024 TYNDP (mean ± min/max)",
    )
    ax.scatter(
        x,
        model_values,
        color="tab:blue",
        zorder=3,
        label="Model",
    )
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Value")
    ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.6)

    if b2a_special:
        levels = ["low", "central", "high"]
        level_colors = {"low": "tab:orange", "central": "tab:green", "high": "tab:red"}
        base_x = len(labels)
        labels.append("B2a_societal_cost_variation")
        ax.set_xticks(list(range(len(labels))))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        offsets = {"low": -0.15, "central": 0.0, "high": 0.15}
        for level in levels:
            model_val = select_value_by_subindex(df, b2a_special, "model", level)
            bench_val = select_value_by_subindex(df, b2a_special, "2024 tyndp", level)
            if bench_val is not None:
                ax.scatter(
                    [base_x + offsets[level]],
                    [bench_val],
                    color=level_colors[level],
                    marker="x",
                    label=f"2024 TYNDP {level}",
                    zorder=3,
                )
            if model_val is not None:
                ax.scatter(
                    [base_x + offsets[level]],
                    [model_val],
                    color=level_colors[level],
                    marker="o",
                    label=f"Model {level}",
                    zorder=4,
                )
    if project_label:
        ax.set_title(f"CBA indicator benchmark ({project_label})")
    else:
        ax.set_title("CBA indicator benchmark")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
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
