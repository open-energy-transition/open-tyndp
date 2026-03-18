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


def create_plots(
    df: pd.DataFrame,
    output_file: str,
    planning_horizons=None,
    area: str = None,
):
    """Create benchmark plot from all collected indicator files."""
    #

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

    plot_items_mEuro = []
    plot_items_percent = []
    axis_items = []
    for planning_horizon in planning_horizons:
        for cyear in cyears:
            for project_id in project_ids:
                # Collect data to plot
                project_label = f"{planning_horizon}, {cyear}, t{project_id}"
                project_df = df[df["project_id"] == project_id].copy()
                indicators = sorted(project_df["indicator"].dropna().unique())

                # Loop through necessary indicators
                for indicator in INDICATOR_UNITS:
                    model_val = select_model_value(project_df, indicator)
                    bench = benchmark_range(
                        project_df, indicator, "TYNDP 2024", planning_horizon
                    )
                    if model_val is None or bench is None:
                        continue
                    plot_items_mEuro.append(bench[1] - model_val)
                    plot_items_percent.append(bench[1] / model_val)
                    axis_items.append(project_label)

    if not plot_items_percent:
        logger.info("No benchmark data to plot")
        return

    plt.figure(figsize=(6.0, 0.3 * len(plot_items_percent)))
    plt.title("CBA differencial cost indicator benchmark\n", y=0.98)
    plt.barh(axis_items, plot_items_percent)
    plt.xlabel("Mean difference B1_total_system_cost_change (%)")
    plt.xticks(rotation=90, fontsize=6)
    plt.yticks(fontsize=6)
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    plt.savefig(output_file, dpi=400)
    plt.close()

    logger.info("Benchmark plots saved to %s", output_file)


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
