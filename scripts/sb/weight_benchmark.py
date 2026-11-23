# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Weights the benchmark data of the weather-cluster runs of a scenario for one planning horizon.

Every row of the benchmark data (Open-TYNDP statistics or TYNDP dashboard outputs) is an
annual quantity of one weather scenario: an energy, a capacity, a cross-border flow, a
time-averaged price or a number of hours. Summing the rows of the selected weather
scenarios with their weather cluster weights, renormalized over the selection, gives the
expected value across the weather scenarios. A gap propagates: a row with a NaN value in,
or missing from, any weather scenario is NaN in the weighted data.

Each price table `{table}` is replaced by two variants with explicit names:

- `{table}_cluster_weighted`: the prices averaged with the weather cluster weights.
- `{table}_load_weighted`: the prices averaged with the weather cluster weights times the
  load of the bus in each weather scenario, taken from the table `{table}_load`.

Outputs
-------
- `results/{run}/benchmarks/tyndp-2026/resources/weighted/{name}{horizon}.csv`: weighted benchmark data in long format.
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)


def weight_benchmarks(
    benchmarks: list[pd.DataFrame], weights: list[float]
) -> pd.DataFrame:
    """
    Combine long-format benchmark data of several weather scenarios with weights.

    Parameters
    ----------
    benchmarks : list[pd.DataFrame]
        Benchmark data per weather scenario with a `value` column.
    weights : list[float]
        Weight of each weather scenario, summing to one.

    Returns
    -------
    pd.DataFrame
        Weighted sum of `value` per combination of the other columns, NaN for a row
        with a NaN value in or missing from any weather scenario. Each price table
        `{table}` is given as `{table}_cluster_weighted` and `{table}_load_weighted`.
    """
    columns = benchmarks[0].columns
    keys = columns.drop("value").tolist()
    values = pd.concat(
        [df.set_index(keys)["value"] for df in benchmarks], axis=1, ignore_index=True
    )
    tables = values.index.get_level_values("table")
    is_load = tables.str.endswith("_load")
    is_price = tables.str.contains("price") & ~is_load
    prices = values[is_price]
    loads = (
        values[is_load]
        .rename(lambda t: t.removesuffix("_load"), level="table")
        .droplevel("unit")
        .reindex(prices.index.droplevel("unit"))
        .set_axis(prices.index)
    )
    weighted = values.dot(weights)
    load_weighted = (prices * loads).dot(weights) / loads.dot(weights)
    return pd.concat(
        [
            weighted[~is_price],
            weighted[is_price].rename(lambda t: f"{t}_cluster_weighted", level="table"),
            load_weighted.rename(lambda t: f"{t}_load_weighted", level="table"),
        ]
    ).reset_index(name="value")[columns]


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "weight_benchmark",
            run="NT",
            name="benchmarks_",
            horizon="2030",
            configfiles="config/config.tyndp.yaml",
        )

    configure_logging(snakemake)

    weights = snakemake.params.weights
    logger.info(f"Weighting {len(weights)} weather scenarios with {weights}")
    benchmarks = [pd.read_csv(fn) for fn in snakemake.input]
    weight_benchmarks(benchmarks, weights).to_csv(snakemake.output[0], index=False)
