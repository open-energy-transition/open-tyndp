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
        benchmark = benchmark[
            benchmark["subindex"].isin(["low", "central", "high"])
        ]
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


def create_plots(df, output_file, area):
    """Create benchmark plot from all collected indicator files."""
    # planning_horizon=None, area=None

    if df.empty:
        logger.warning("No indicators data to plot")
        return

    output_path = Path(output_file).parent
    indicators = df.indicator.unique()
    subindexes = df.subindex.unique()
    planning_horizons = df.planning_horizon.unique()
    project_ids = df.loc[df["source"] == "Open-TYNDP", "project_id"].dropna().unique()
    cyears = df.loc[df["source"] == "Open-TYNDP", "cyear"].dropna().unique()

    if len(planning_horizons) == 0:
        logger.info("No planning horizons found in dataset")
        return
    if len(project_ids) == 0:
        logger.info("No model projects found in dataset")
        return
    if len(indicators) == 0:
        logger.info("No benchmark indicators available in dataset")
        return
    if len(subindexes) == 0:
        logger.info("No benchmark subindexes available in dataset")
        return

    plot_items = []
    axis_items = []
    for planning_horizon in planning_horizons:
        for cyear in cyears:
            for project_id in project_ids:
                # Collect data to plot
                project_label = f'{planning_horizon}_{cyear}_{project_id}'
                project_df = df[df["project_id"] == project_id].copy()
                indicators = sorted(project_df["indicator"].dropna().unique())

                # Loop through necessary indicators
                for indicator in INDICATOR_UNITS:
                    model_val = select_model_value(project_df, indicator)
                    bench = benchmark_range(project_df, indicator, "TYNDP 2024", planning_horizon)
                    if model_val is None or bench is None:
                        continue
                    plot_items.append(bench[1] - model_val)
                    axis_items.append(project_label)

    if not plot_items:
        logger.info("No benchmark data to plot")
        return

    plt.figure(figsize=(10, 4.2))
    plt.bar(axis_items, plot_items)
    plt.xticks(rotation=90)
    plt.savefig(output_file, dpi=400)
    plt.close()

    """
    fig, axes = plt.subplots(
        nrows=1,
        ncols=1,
        figsize=(max(1.7 * len(plot_items), 10), 4.2),
        squeeze=False,
        sharey=False,
    )
    axes[0].bar(axis_items, plot_items)
    axes[0].set_ylabel('B1 delta (mEuro/year)')
    axes[0].set_title('B1 Cost Delta per run')
    fig.tight_layout(rect=[0, 0.12, 1, 0.90])
    fig.savefig(output_path, dpi=400)
    plt.close(fig)
    """

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
    
    # Create a list of files
    input_files = str(snakemake.input).split(" ")
    output_file = str(snakemake.output)

    # Read all files into one DataFrame
    df = pd.concat(map(pd.read_csv, input_files), ignore_index=True)
    area = snakemake.config.get("cba", {}).get("area")

    # Create a summary plot for a specific project
    create_plots(df, output_file, area)
