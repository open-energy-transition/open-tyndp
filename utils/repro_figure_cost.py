# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Time each stage of one benchmark figure to locate a per-figure slowdown.

Mirrors the figure built by `_plot_scenario_comparison` in
`scripts/sb/plot_benchmark.py`, so a stage that dominates here is the stage
that dominates a real `plot_benchmark` run. Temporary diagnostic; run from the
repository root:

    pixi run -e open-tyndp python utils/repro_figure_cost.py
"""

import os
import sys
import time
from collections.abc import Iterator
from contextlib import contextmanager

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties, findfont

plt.style.use("bmh")

TIMINGS: list[tuple[str, float]] = []

N_CARRIERS = 14  # worst real figure, measured from benchmarks_s_all___*.csv
SOURCES = [
    "Open-TYNDP",
    "TYNDP 2024 Scenarios Report",
    "TYNDP 2024 Vis Pltfm",
    "Market Model",
]


@contextmanager
def timed(name: str) -> Iterator[None]:
    """Record the wall-clock duration of the enclosed block under `name`."""
    start = time.perf_counter()
    yield
    TIMINGS.append((name, time.perf_counter() - start))


def report_environment() -> None:
    """Print the versions and paths that determine rendering cost."""
    print("--- environment ---")
    print(f"python      : {sys.version.split()[0]}")
    print(f"matplotlib  : {matplotlib.__version__}  backend={matplotlib.get_backend()}")
    print(f"cache dir   : {matplotlib.get_cachedir()}")
    print(f"numpy       : {np.__version__}   pandas: {pd.__version__}")
    try:
        from matplotlib.ft2font import __freetype_version__

        print(f"freetype    : {__freetype_version__}")
    except ImportError as error:
        print(f"freetype    : unavailable ({error})")


def sample_frame() -> pd.DataFrame:
    """Build a DataFrame shaped like the widest real benchmark figure."""
    values = np.random.default_rng(0).random((N_CARRIERS, len(SOURCES))) * 500
    return pd.DataFrame(
        values,
        index=[f"carrier {i} with a fairly long name" for i in range(N_CARRIERS)],
        columns=SOURCES,
    )


def time_font_lookup() -> None:
    """Time a first font resolution, which populates the font cache."""
    with timed("findfont (first)"):
        font = findfont(FontProperties(family=plt.rcParams["font.family"]))
    print(f"font.family : {plt.rcParams['font.family']} -> {os.path.basename(font)}")


def time_get_version() -> None:
    """Time `get_version`, which `add_metadata` calls once per figure."""
    sys.path.insert(0, os.getcwd())
    try:
        from scripts._helpers import get_version
    except ImportError as error:
        print(f"get_version : unavailable ({error})")
        return
    with timed("get_version (git)"):
        get_version()


def time_one_figure(df: pd.DataFrame) -> None:
    """Time every stage of building, laying out and writing one figure."""
    with timed("Figure + canvas"):
        fig = Figure(figsize=(12, 8))
        FigureCanvasAgg(fig)
        ax = fig.subplots()

    with timed("plot.bar"):
        df.clip(lower=0).plot.bar(
            ax=ax, width=0.7, xlabel="", ylabel="Power Capacity [GW]", title="bench"
        )

    with timed("tick labels"):
        ax.tick_params(axis="x", labelrotation=45)
        plt.setp(ax.get_xticklabels(), ha="right")

    with timed("legend (explicit loc)"):
        ax.legend(loc="upper left", facecolor="white", frameon=True, edgecolor="black")

    with timed(f"bar_label x{len(ax.containers)}"):
        for container in ax.containers:
            ax.bar_label(container, fmt=lambda x: f"{x:.1f}", padding=3, fontsize=8)

    time_get_version()

    with timed("draw_without_rendering"):
        fig.draw_without_rendering()

    with timed("get_tightbbox"):
        fig.get_tightbbox(fig.canvas.get_renderer())

    with timed("savefig pdf (plain)"):
        fig.savefig(os.devnull, format="pdf")

    with timed("savefig pdf (tight)"):
        fig.savefig(os.devnull, format="pdf", bbox_inches="tight")

    with timed("savefig png"):
        fig.savefig(os.devnull, format="png", dpi=100)


def report_timings() -> None:
    """Print recorded stage durations, slowest first."""
    total = sum(duration for _, duration in TIMINGS)
    print()
    for name, duration in sorted(TIMINGS, key=lambda item: -item[1]):
        share = duration / total * 100 if total else 0.0
        print(f"  {name:26s} {duration * 1000:9.1f} ms   {share:5.1f}%")
    print(f"  {'TOTAL':26s} {total * 1000:9.1f} ms")


def main() -> None:
    """Report the environment, then time one figure stage by stage."""
    report_environment()
    time_font_lookup()
    print("\n--- one figure, staged (mirrors _plot_scenario_comparison) ---")
    time_one_figure(sample_frame())
    report_timings()


if __name__ == "__main__":
    main()
