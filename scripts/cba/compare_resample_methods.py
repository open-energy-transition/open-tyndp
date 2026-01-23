# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Compare MSV resampling methods for CBA networks.

This script analyzes the impact of different resampling methods (ffill, bfill,
nearest, interpolate) when upsampling Marginal Storage Values from coarse to
fine temporal resolution.

The analysis outputs:
- objective_comparison.csv: System cost differences between methods
- msv_statistics.csv: MSV statistics per method and store
- msv_pairwise_differences.csv: Max/mean absolute differences vs reference
- msv_snapshot_differences.csv: Count of differing snapshots per store
"""

import logging
from pathlib import Path

import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def load_networks(network_paths: dict[str, str]) -> dict[str, pypsa.Network]:
    """Load networks from paths dictionary."""
    networks = {}
    for name, path in network_paths.items():
        logger.info(f"Loading {name} network from {path}")
        networks[name] = pypsa.Network(path)
    return networks


def extract_msv_comparison(networks: dict[str, pypsa.Network]) -> pd.DataFrame:
    """Extract time-varying MSV from all networks into a comparison DataFrame."""
    dfs = []
    for method, n in networks.items():
        msv = n.stores_t.marginal_cost.copy()
        msv = msv.loc[:, (msv != 0).any()]
        msv_melted = msv.reset_index().melt(
            id_vars="snapshot", var_name="store", value_name="msv"
        )
        msv_melted["method"] = method
        dfs.append(msv_melted)
    return pd.concat(dfs, ignore_index=True)


def compute_objective_summary(networks: dict[str, pypsa.Network]) -> pd.DataFrame:
    """Compute objective function comparison across methods."""
    data = {
        method: {"objective_eur": n.objective, "snapshots": len(n.snapshots)}
        for method, n in networks.items()
    }
    df = pd.DataFrame(data).T
    df.index.name = "method"

    base = df.loc["ffill", "objective_eur"]
    df["diff_from_ffill_eur"] = df["objective_eur"] - base
    df["diff_from_ffill_pct"] = (df["objective_eur"] - base) / base * 100
    return df.reset_index()


def compute_msv_statistics(msv_df: pd.DataFrame) -> pd.DataFrame:
    """Compute MSV statistics per method and store."""
    return (
        msv_df.groupby(["method", "store"])["msv"]
        .agg(["mean", "std", "min", "max"])
        .reset_index()
    )


def compute_pairwise_differences(
    networks: dict[str, pypsa.Network], reference: str = "ffill"
) -> pd.DataFrame:
    """Compute max absolute MSV differences between methods."""
    ref_msv = networks[reference].stores_t.marginal_cost
    stores = ref_msv.columns[ref_msv.abs().sum() > 0]

    records = []
    for method, n in networks.items():
        if method == reference:
            continue
        msv = n.stores_t.marginal_cost[stores]
        diff = (ref_msv[stores] - msv).abs()
        for store in stores:
            records.append(
                {
                    "reference": reference,
                    "method": method,
                    "store": store,
                    "max_abs_diff": diff[store].max(),
                    "mean_abs_diff": diff[store].mean(),
                }
            )
    return pd.DataFrame(records)


def compute_snapshot_differences(
    networks: dict[str, pypsa.Network], reference: str = "ffill"
) -> pd.DataFrame:
    """Count snapshots with differences per store."""
    ref_msv = networks[reference].stores_t.marginal_cost
    stores = ref_msv.columns[ref_msv.abs().sum() > 0]

    records = []
    for method, n in networks.items():
        if method == reference:
            continue
        msv = n.stores_t.marginal_cost[stores]
        for store in stores:
            n_diff = (ref_msv[store] != msv[store]).sum()
            records.append(
                {
                    "reference": reference,
                    "method": method,
                    "store": store,
                    "snapshots_different": n_diff,
                    "total_snapshots": len(ref_msv),
                    "pct_different": n_diff / len(ref_msv) * 100,
                }
            )
    return pd.DataFrame(records)


def run_comparison(network_paths: dict[str, str], output_dir: str):
    """Run comparison analysis and save results."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    networks = load_networks(network_paths)

    logger.info("Computing objective comparison")
    compute_objective_summary(networks).to_csv(
        output_path / "objective_comparison.csv", index=False
    )

    logger.info("Computing MSV statistics")
    msv_df = extract_msv_comparison(networks)
    compute_msv_statistics(msv_df).to_csv(
        output_path / "msv_statistics.csv", index=False
    )

    logger.info("Computing pairwise differences")
    compute_pairwise_differences(networks).to_csv(
        output_path / "msv_pairwise_differences.csv", index=False
    )

    logger.info("Computing snapshot differences")
    compute_snapshot_differences(networks).to_csv(
        output_path / "msv_snapshot_differences.csv", index=False
    )

    logger.info(f"Results saved to {output_path}")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("compare_msv_resample_methods")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    network_paths = {
        "ffill": snakemake.input.ffill,
        "nearest": snakemake.input.nearest,
        "interpolate": snakemake.input.interpolate,
        "bfill": snakemake.input.bfill,
    }

    run_comparison(network_paths, snakemake.output.output_dir)
