# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Collect indicators from individual indicators files and create weighted average
indicator.

This script reads all collected project indicators from the different weather years
and calculates a new combined indicator which is stores in a new CSV file.
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

INDICATOR_UNITS = {
    "B1_total_system_cost_change": "Meuro/year",
    #    "B2a_co2_variation": "t/year",
    #    "B2a_societal_cost_variation": "Meuro/year",
    #    "B3a_res_capacity_change": "MW",
    #    "B3_res_generation_change": "MWh/year",
    #    "B3_annual_avoided_curtailment": "MWh/year",
    #    "B4a_nox": "kg/year",
    #    "B4b_nh3": "kg/year",
    #    "B4c_sox": "kg/year",
    #    "B4d_pm25": "kg/year",
    #    "B4e_pm10": "kg/year",
    #    "B4f_nmvoc": "kg/year",
}

# TODO read from CSV file
cyear_weightings = {
    1995: 0.233,
    2008: 0.367,
    2009: 0.400,
}


def select_model_value(df: pd.DataFrame, indicator: str) -> float | None:
    """Pick a representative model value for a given indicator."""
    model = df[(df["source"] == "Open-TYNDP") & (df["indicator"] == indicator)]
    if model.empty:
        return None

    for subindex in ["mean", "central", ""]:
        subset = model[model["subindex"] == subindex]
        if not subset.empty:
            return subset["value"].mean()

    return model["value"].mean()


def select_value_by_subindex(
    df: pd.DataFrame, indicator: str, source: str, subindex: str, planning_horizon: int
) -> float | None:
    """Return the mean value for a given indicator/subindex/source."""
    subset = df[
        (df["source"] == source)
        & (df["indicator"] == indicator)
        & (df["subindex"] == subindex)
        & (df["planning_horizon"] == str(planning_horizon))
    ]
    if subset.empty:
        return None
    return float(subset["value"].mean())


def benchmark_range(
    df: pd.DataFrame,
    indicator: str,
    source: str = "TYNDP 2024",
    planning_horizon: int = 0,
) -> tuple[float, float, float] | None:
    """Return (min, mean, max) range for a benchmark indicator."""
    benchmark = df[
        (df["source"] == source)
        & (df["indicator"] == indicator)
        & (df["planning_horizon"] == planning_horizon)
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


def format_area_subtitle(area: str | None) -> str | None:
    """Map cba.area to a plot subtitle."""
    if not area:
        return None
    mapping = {
        "tyndp": "Impact on whole area of TYNDP Study",
        "entso-e": "Impact on ENTSO-E area",
        "eu27": "Impact on EU27 area",
    }
    return mapping.get(str(area).lower())


def plot_project_benchmarks(
    df: pd.DataFrame,
    output_path: Path,
    project_label: str | None = None,
    area_subtitle: str | None = None,
) -> None:
    """Plot one subplot per indicator with its own y-axis and legend."""
    indicators = sorted(df["indicator"].dropna().unique())
    if not indicators:
        logger.info("No benchmark indicators available to plot")
        return

    planning_horizons = sorted(df["planning_horizon"].dropna().unique())
    if not planning_horizons:
        logger.info("No benchmark planning horizon available to plot")
        return

    plot_items = []
    # loop through necessary indicators
    for indicator in INDICATOR_UNITS:
        model_val = select_model_value(df, indicator)
        bench = benchmark_range(df, indicator, "TYNDP 2024", planning_horizons[0])
        if model_val is None or bench is None:
            continue
        plot_items.append(indicator)

    if not plot_items:
        logger.info("No complete benchmark data to plot")
        return

    fig, axes = plt.subplots(
        nrows=1,
        ncols=len(plot_items) * len(planning_horizons),
        figsize=(max(1.7 * len(plot_items), 6), 4.2),
        squeeze=False,
        sharey=True,
    )

    legend_handles = []
    legend_labels = []
    for ax, planning_horizon in zip(axes[0], planning_horizons):
        ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.6)
        ax.set_title(planning_horizon, fontstyle="italic")
        model_range = benchmark_range(df, indicator, "Open-TYNDP", planning_horizon)
        if model_range is None:
            model_min_val, model_mean_val, model_max_val = (0, 0, 0)
        else:
            model_min_val, model_mean_val, model_max_val = model_range
        min_val, mean_val, max_val = benchmark_range(
            df, indicator, "TYNDP 2024", planning_horizon
        )
        ax.errorbar(
            [-0.1],
            [mean_val],
            yerr=[[abs(mean_val - min_val)], [abs(max_val - mean_val)]],
            fmt="x",
            color="gray",
            ecolor="lightgray",
            capsize=3,
        )
        label = "2024 TYNDP (mean ± min/max)"
        if label not in legend_labels:
            legend_handles.append(
                Line2D([0], [0], marker="x", color="gray", linestyle="None")
            )
            legend_labels.append(label)
        ax.errorbar(
            [0.1],
            [model_mean_val],
            yerr=[
                [abs(model_mean_val - model_min_val)],
                [abs(model_max_val - model_mean_val)],
            ],
            fmt="o",
            color="tab:blue",
            ecolor="lightblue",
            capsize=3,
        )
        ax.set_xlim(xmin=-0.5, xmax=0.5)
        label = "Open-TYNDP (mean ± min/max)"
        if label not in legend_labels:
            legend_handles.append(
                Line2D([0], [0], marker="o", color="tab:blue", linestyle="None")
            )
            legend_labels.append(label)
        ax.set_xticks([])
        ylim = ax.get_ylim()
        values = [v for v in ylim if v != 0]
        if values:
            abs_max = max(abs(v) for v in values)
            abs_min = min(abs(v) for v in values)
            if abs_min > 0 and abs_max / abs_min >= 1e3:
                ax.set_yscale("symlog", linthresh=max(abs_min, 1.0))

        units = df.loc[
            (df["indicator"] == indicator) & (df["source"] == "Open-TYNDP"), "units"
        ].dropna()
        if planning_horizons.index(planning_horizon) == 0:
            unit_label = units.iloc[0] if not units.empty else ""
            title = f"{indicator} ({unit_label})" if unit_label else indicator
            ax.set_ylabel(title)

    if project_label:
        fig.suptitle(f"CBA indicator benchmark ({project_label})", y=0.98)
    if area_subtitle:
        fig.text(0.5, 0.92, area_subtitle, ha="center", va="center")
    # Interleaved to account for matplotlib's column-major legend layout.
    other_items = [(h, l) for h, l in zip(legend_handles, legend_labels)]

    if other_items:
        other_handles = [h for h, _ in other_items]
        other_labels = [l for _, l in other_items]
        fig.legend(
            other_handles,
            other_labels,
            loc="lower center",
            bbox_to_anchor=(0.42, -0.02),
            ncol=1,
            frameon=False,
        )
    fig.tight_layout(rect=[0, 0.12, 1, 0.90])
    fig.savefig(output_path, dpi=400)
    plt.close(fig)


def create_plots(df, output_file, project_id, area):
    """Create benchmark plots from a per-project or collected indicators file."""
    # planning_horizon=None, area=None

    if df.empty:
        logger.warning("No indicators data to plot")
        return

    planning_horizon = df.planning_horizon.unique()[0]
    output_path = Path(output_file).parent
    project_ids = df.loc[df["source"] == "Open-TYNDP", "project_id"].dropna().unique()
    if output_path.suffix:
        if len(project_ids) == 0:
            logger.info("No model projects found in indicators file")
            return
        project_id = int(project_ids[0])
        project_code_series = df.loc[
            df["project_id"] == project_id, "project_code"
        ].dropna()
        project_code = str(project_code_series.iloc[0])
        project_label = project_code
        project_df = df[df["project_id"] == project_id].copy()
        plot_project_benchmarks(
            project_df, output_path, project_label, format_area_subtitle(area)
        )
    else:
        for project_id in project_ids:
            project_df = df[df["project_id"] == project_id].copy()
            project_id = int(project_id)
            project_code_series = df.loc[
                df["project_id"] == project_id, "project_code"
            ].dropna()
            project_code = str(project_code_series.iloc[0])
            suffix = f"_{planning_horizon}" if planning_horizon else ""
            project_label = project_code
            plot_project_benchmarks(
                project_df, output_file, project_label, format_area_subtitle(area)
            )

    logger.info("Benchmark plots saved to %s", output_file)


def summarize_indicators(input_files, output_file):
    """
    Concatenate multiple CSV files into one using the csv module.

    Args:
        input_files: List of paths to input CSV files
        output_file: Path to output CSV file

    The function:
    1. Reads the header from the first file
    2. Writes all rows from all files to the output
    3. Ensures all files have the same header structure
    """
    if not input_files:
        logger.warning("No input files provided")
        # Create empty output file with no header
        with open(output_file, "w", newline="") as f:
            pass
        return

    # Read input files
    logger.info(f"Read {len(input_files)} indicator files")
    df = pd.concat(map(pd.read_csv, ["mydata.csv", "mydata1.csv"]), ignore_index=True)

    project_ids = df.project_id.unique()[0]
    alternative_subindex = {"min": "low", "mean": "central", "max": "high"}

    # loop through all coolected indicators
    for INDICATOR_UNIT in INDICATOR_UNITS:
        units = df[(df.indicator == INDICATOR_UNIT)].units.unique()[0]
        indicator_values = {
            "min": (
                df[(df.indicator == INDICATOR_UNIT) & (df.source == "Open-TYNDP")].value
            ).min(),
            "mean": (
                df[(df.indicator == INDICATOR_UNIT) & (df.source == "Open-TYNDP")].value
                * df[
                    (df.indicator == INDICATOR_UNIT) & (df.source == "Open-TYNDP")
                ].cyear_weight
            ).sum(),
            "max": (
                df[(df.indicator == INDICATOR_UNIT) & (df.source == "Open-TYNDP")].value
            ).max(),
        }


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("collect_indicators")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Validate input is available
    if not snakemake.input:
        raise ValueError("Reference network is not solved")

    # Read input files
    logger.info(f"Read {len(snakemake.input)} indicator files")

    # create a list of files
    input_files = str(snakemake.input).split(" ")
    output_file = str(snakemake.output)
    # read all files into one DataFrame
    df = pd.concat(map(pd.read_csv, input_files), ignore_index=True)
    project_id = df.project_id.unique()[0]
    area = snakemake.config.get("cba", {}).get("area")

    # Create a summary plot for a specific project
    create_plots(df, output_file, project_id, area)
