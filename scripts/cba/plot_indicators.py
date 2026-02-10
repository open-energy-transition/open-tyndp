# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Create plots for CBA indicators.

This script reads the collected indicators CSV file and generates various
plots to visualize the cost-benefit analysis results, including the B1
indicator (Total System Cost difference).
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def save_figure(fig, output_dir, filename, output_formats):
    """Save figure in specified format(s)."""
    for fmt in output_formats:
        filepath = output_dir / f"{filename}.{fmt}"
        save_kwargs = {"bbox_inches": "tight"}
        if fmt == "png":
            save_kwargs["dpi"] = 150
        fig.savefig(filepath, **save_kwargs)
    logger.info(f"Saved {filename}")


def get_indicator_series(df: pd.DataFrame, indicator: str) -> pd.Series:
    """Return a project_id-indexed Series for a single indicator."""
    series = df.loc[df["indicator"] == indicator, ["project_id", "value"]]
    if series.empty:
        return pd.Series(dtype=float)
    return series.set_index("project_id")["value"]


def load_and_merge_data(indicators_path, projects_path):
    """Load indicators and merge with project metadata."""
    raw = pd.read_csv(indicators_path)
    projects = pd.read_csv(projects_path)

    border_counts = projects.groupby("project_id").size().rename("border_count")
    projects_agg = (
        projects.groupby("project_id")
        .agg(
            {
                "project_name": "first",
                "is_crossborder": "first",
                "p_nom 0->1": "sum",
                "p_nom 1->0": "sum",
            }
        )
        .reset_index()
    )
    projects_agg = projects_agg.merge(border_counts, on="project_id")
    projects_agg["total_capacity_MW"] = (
        projects_agg["p_nom 0->1"] + projects_agg["p_nom 1->0"]
    )

    base = raw.copy()
    # filter for Open-TYNDP source (to get model results)
    if "source" in base.columns:
        base = base[base["source"] == "Open-TYNDP"]

    meta_cols = [c for c in ["method", "is_beneficial", "interpretation"] if c in base]
    meta = base.groupby("project_id")[meta_cols].first().reset_index()
    print(meta)

    metrics = meta[["project_id"]].copy()
    metrics["B1_total_system_cost_change"] = metrics["project_id"].map(
        get_indicator_series(base, "B1_total_system_cost_change")
    )
    metrics["capex_change"] = metrics["project_id"].map(
        get_indicator_series(base, "capex_change")
    )
    metrics["opex_change"] = metrics["project_id"].map(
        get_indicator_series(base, "opex_change")
    )
    indicators_wide = meta.merge(metrics, on="project_id", how="left")

    merged = indicators_wide.merge(projects_agg, on="project_id", how="left")
    if "is_beneficial" in merged.columns:
        merged["is_beneficial"] = (
            merged["is_beneficial"]
            .astype(str)
            .str.lower()
            .map({"true": True, "false": False, "1": True, "0": False})
        )
    merged["B1_billion_EUR"] = merged["B1_total_system_cost_change"] / 1e9
    merged["capex_change_billion"] = merged["capex_change"] / 1e9
    merged["opex_change_billion"] = merged["opex_change"] / 1e9

    return merged, projects["project_id"].nunique()


def plot_b1_top_projects(df, output_dir, method, colors, output_formats, n_top=20):
    """
    Diverging bar chart showing B1 with CAPEX/OPEX breakdown for top N projects.

    Each project is shown on a single row with:
    - CAPEX change extending left from zero (negative = cost saved)
    - OPEX change extending right from zero (positive = additional cost)
    - Diamond marker showing the net B1 value
    """
    df_top = df.nlargest(n_top, "B1_billion_EUR", keep="first")
    df_sorted = df_top.sort_values("B1_billion_EUR", ascending=True).reset_index(
        drop=True
    )

    fig, ax = plt.subplots(figsize=(14, max(6, n_top * 0.5)))
    y_pos = np.arange(len(df_sorted))
    bar_height = 0.35

    ax.barh(
        y_pos,
        df_sorted["capex_change_billion"],
        height=bar_height,
        color=colors["capex"],
        alpha=0.9,
        edgecolor="white",
        linewidth=0.5,
    )
    ax.barh(
        y_pos,
        df_sorted["opex_change_billion"],
        height=bar_height,
        color=colors["opex"],
        alpha=0.9,
        edgecolor="white",
        linewidth=0.5,
    )

    for i, (_, row) in enumerate(df_sorted.iterrows()):
        color = (
            colors["beneficial"]
            if row["B1_billion_EUR"] >= 0
            else colors["not_beneficial"]
        )
        ax.plot(
            row["B1_billion_EUR"],
            i,
            "D",
            color=color,
            markersize=10,
            markeredgecolor="white",
            markeredgewidth=2,
            zorder=5,
        )

    labels = []
    for _, row in df_sorted.iterrows():
        multi = "*" if row["border_count"] > 1 else ""
        labels.append(
            f"{row['project_id']}: {row['project_name']}{multi} "
            f"[{row['p_nom 0->1']:.0f}/{row['p_nom 1->0']:.0f} MW]"
        )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=10)

    ax.axvline(x=0, color="black", linewidth=1.2, zorder=3)
    ax.set_xlabel(
        "◄ Savings | Costs ►  Cost Change (Billion EUR)",
        fontsize=12,
        fontweight="bold",
    )
    ax.set_title(
        f"{method.upper()} B1 Indicator: Top {n_top} Projects by Benefit",
        fontsize=14,
        fontweight="bold",
    )

    max_val = max(
        df_sorted["B1_billion_EUR"].abs().max(),
        df_sorted["opex_change_billion"].abs().max(),
        df_sorted["capex_change_billion"].abs().max(),
    )
    for i, (_, row) in enumerate(df_sorted.iterrows()):
        total = row["B1_billion_EUR"]
        color = colors["beneficial"] if total >= 0 else colors["not_beneficial"]
        ax.text(
            max_val * 1.15,
            i,
            f"{total:+.2f}B",
            va="center",
            ha="left",
            fontsize=9,
            color=color,
            fontweight="medium",
        )

    ax.grid(axis="x", alpha=0.3, linestyle="--")
    ax.set_axisbelow(True)

    legend_elements = [
        plt.Rectangle(
            (0, 0), 1, 1, facecolor=colors["capex"], alpha=0.9, label="CAPEX Change"
        ),
        plt.Rectangle(
            (0, 0), 1, 1, facecolor=colors["opex"], alpha=0.9, label="OPEX Change"
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            color="w",
            markerfacecolor=colors["beneficial"],
            markersize=10,
            markeredgecolor="white",
            markeredgewidth=2,
            label="B1 (Beneficial)",
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            color="w",
            markerfacecolor=colors["not_beneficial"],
            markersize=10,
            markeredgecolor="white",
            markeredgewidth=2,
            label="B1 (Not Beneficial)",
        ),
    ]
    ax.legend(
        handles=legend_elements,
        loc="upper left",
        bbox_to_anchor=(0, -0.12),
        ncol=4,
        fontsize=9,
        frameon=False,
    )
    ax.set_xlim(ax.get_xlim()[0], max_val * 1.4)

    plt.tight_layout()
    save_figure(fig, output_dir, f"b1_top_{n_top}_projects", output_formats)
    plt.close(fig)


def plot_b1_summary(df, output_dir, method, colors, output_formats, total_projects):
    """Summary plot with B1 histogram."""
    beneficial = df[df["is_beneficial"] == True]
    not_beneficial = df[df["is_beneficial"] == False]

    fig, ax = plt.subplots(figsize=(10, 6))

    bins = np.linspace(
        df["B1_billion_EUR"].min() - 0.1, df["B1_billion_EUR"].max() + 0.1, 30
    )
    ax.hist(
        beneficial["B1_billion_EUR"],
        bins=bins,
        color=colors["beneficial"],
        alpha=0.7,
        label=f"Beneficial ({len(beneficial)})",
        edgecolor="white",
    )
    ax.hist(
        not_beneficial["B1_billion_EUR"],
        bins=bins,
        color=colors["not_beneficial"],
        alpha=0.7,
        label=f"Not Beneficial ({len(not_beneficial)})",
        edgecolor="white",
    )

    ax.axvline(x=0, color="black", linewidth=1.5, linestyle="--")
    ax.set_xlabel("B1 Indicator (Billion EUR)", fontsize=11)
    ax.set_ylabel("Number of Projects", fontsize=11)
    ax.set_title(
        f"{method.upper()} B1: Distribution of Project Benefits",
        fontsize=12,
        fontweight="bold",
    )
    ax.legend(fontsize=10)
    ax.grid(axis="y", alpha=0.3)

    # Force integer y-axis ticks (number of projects must be integers)
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    plt.tight_layout()
    save_figure(fig, output_dir, "b1_summary", output_formats)
    plt.close(fig)


def plot_b1_capex_vs_opex(df, output_dir, method, colors, output_formats):
    """Scatter plot of B1 CAPEX vs OPEX changes."""
    fig, ax = plt.subplots(figsize=(10, 8))

    beneficial = df[df["is_beneficial"]]
    not_beneficial = df[~df["is_beneficial"]]

    ax.scatter(
        beneficial["capex_change_billion"],
        beneficial["opex_change_billion"],
        c=colors["beneficial"],
        s=80,
        alpha=0.7,
        label=f"Beneficial ({len(beneficial)})",
        edgecolors="white",
        linewidth=0.5,
    )
    ax.scatter(
        not_beneficial["capex_change_billion"],
        not_beneficial["opex_change_billion"],
        c=colors["not_beneficial"],
        s=80,
        alpha=0.7,
        label=f"Not Beneficial ({len(not_beneficial)})",
        edgecolors="white",
        linewidth=0.5,
    )

    ax.axhline(y=0, color="gray", linewidth=0.8, linestyle="--")
    ax.axvline(x=0, color="gray", linewidth=0.8, linestyle="--")

    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1]),
    ]
    ax.plot(lims, [-x for x in lims], "k:", alpha=0.5, label="CAPEX + OPEX = 0")

    ax.set_xlabel("CAPEX Change (Billion EUR)", fontsize=11)
    ax.set_ylabel("OPEX Change (Billion EUR)", fontsize=11)
    ax.set_title(
        f"{method.upper()} B1: CAPEX vs OPEX Changes\n"
        "(Projects above diagonal are beneficial)",
        fontsize=12,
        fontweight="bold",
    )
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    save_figure(fig, output_dir, "b1_capex_vs_opex", output_formats)
    plt.close(fig)


# Future indicator plot functions:
# def plot_b2_...
# def plot_b3_...


def create_plots(indicators_path, projects_path, output_dir, params):
    """Create all CBA indicator plots."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cba_config = params.plotting["cba"]
    colors = cba_config["colors"]
    output_formats = cba_config.get("output_format", ["svg"])
    if isinstance(output_formats, str):
        output_formats = [output_formats]

    df, total_projects = load_and_merge_data(indicators_path, projects_path)
    if df.empty:
        logger.warning("No indicators data to plot")
        return

    method = df["method"].iloc[0].lower()
    logger.info(f"Creating {method.upper()} plots for {len(df)} projects")

    # B1 indicator plots
    plot_b1_top_projects(
        df, output_dir, method, colors, output_formats, n_top=min(20, len(df))
    )
    plot_b1_summary(df, output_dir, method, colors, output_formats, total_projects)
    plot_b1_capex_vs_opex(df, output_dir, method, colors, output_formats)

    # Future indicators (add when implemented in make_indicators.py):
    # plot_b2_...(df, output_dir, method)
    # plot_b3_...(df, output_dir, method)

    logger.info(f"All plots saved to {output_dir}")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_indicators",
            cba_method="toot",
            planning_horizons="2030",
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    create_plots(
        snakemake.input.indicators,
        snakemake.input.transmission_projects,
        snakemake.output.plot_dir,
        snakemake.params,
    )
