# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script computes accuracy indicators for comparing workflow results against reference data from the TYNDP 2024.

This module implements a methodology introduced by Wen et al. (2022) for evaluating performance of energy system
models using multiple accuracy indicators.
"""

import logging
import multiprocessing as mp
import os
from functools import partial

import numpy as np
import pandas as pd
from tqdm import tqdm

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def load_data(benchmarks_fn: str, results_fn: str, scenario: str) -> pd.DataFrame:
    """
    Load Open-TYNDP and TYNDP 2024 results.

    Parameters
    ----------
    benchmarks_fn : str
        Path to the TYNDP 2024 benchmark data file.
    results_fn : str
        Path to the Open-TYNDP results data file.
    scenario: str
        Name of scenario to compare.

    Returns
    -------
    benchmarks_raw : pd.DataFrame
        Combined DataFrame containing both Open-TYNDP and TYNDP 2024 data.

    """

    # Load data
    logger.info("Loading benchmark using TYNDP 2024 and Open-TYNDP")
    benchmarks_tyndp = pd.read_csv(benchmarks_fn).query("scenario==@scenario")
    benchmarks_n = []
    for fn in results_fn:
        benchmarks_n.append(pd.read_csv(fn).query("scenario==@scenario"))
    benchmarks_n = pd.concat(benchmarks_n)
    benchmarks_raw = pd.concat([benchmarks_tyndp, benchmarks_n]).dropna(
        how="all", axis=1
    )

    # Clean data
    available_years = set(benchmarks_tyndp.year).intersection(benchmarks_n.year)  # noqa: F841
    benchmarks_raw = benchmarks_raw.query("year in @available_years")

    return benchmarks_raw


def match_temporal_resolution(
    df: pd.DataFrame,
    model_col: str = "Open-TYNDP",
    rfc_col: str = "TYNDP 2024",
) -> pd.DataFrame:
    """
    Match temporal resolution against reference data. Hourly time series from the rfc_col will be
    aggregated to match the temporal resolution of model_col.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with MultiIndex containing snapshot timestamps and data columns.
    model_col : str, default "Open-TYNDP"
        Column name for model values with potentially lower temporal resolution.
    rfc_col : str, default "TYNDP 2024"
        Column name for reference values with hourly temporal resolution.

    Returns
    -------
    pd.DataFrame
        Aggregated DataFrame where reference data temporal resolution matches model data.
    """
    idx_agg = (
        df[model_col]
        .dropna()
        .index.get_level_values("snapshot")
        .drop_duplicates()
        .sort_values()
    )
    idx_full = (
        df[rfc_col]
        .dropna()
        .index.get_level_values("snapshot")
        .drop_duplicates()
        .sort_values()
    )

    aggregation_map = (
        pd.Series(idx_agg.rename("map"), index=idx_agg).reindex(idx_full).ffill()
    )
    df_map = df.join(aggregation_map, on="snapshot", how="left")

    df_agg = (
        df_map.groupby(["carrier", "scenario", "year", "table", "map"])
        .mean()
        .rename_axis(index={"map": "snapshot"})
    )

    return df_agg


def _compute_smpe(df: pd.DataFrame, model_col: str, rfc_col: str, eps: float) -> float:
    """
    Calculate Symmetric Mean Percentage Error (sMPE).

    sMPE indicates the direction of the deviations between modeled scenarios
    and reference outcomes, showing if the output is overall overestimated or underestimated.

    Formula: sMPE = (1/n) * Σ[(ŷᵢ - yᵢ) / ((|ŷᵢ| + |yᵢ|)/2 + ε)]

    Reference
    ---------
    Wen et al. (2022), Applied Energy 325, 119906, Table 1
    """
    return (
        (df[model_col] - df[rfc_col])
        / ((df[model_col].abs() + df[rfc_col].abs()) / 2 + eps)
    ).mean()


def _compute_smape(df: pd.DataFrame, model_col: str, rfc_col: str, eps: float) -> float:
    """
    Calculate Symmetric Mean Absolute Percentage Error (sMAPE).

    sMAPE indicates the absolute magnitude of the deviations. It avoids the cancellation of negative and positive errors.

    Formula: sMAPE = (1/n) * Σ[|ŷᵢ - yᵢ| / ((|ŷᵢ| + |yᵢ|)/2 + ε)]

    Reference
    ---------
    Wen et al. (2022), Applied Energy 325, 119906, Table 1
    """
    return (
        (df[model_col] - df[rfc_col]).abs()
        / ((df[model_col].abs() + df[rfc_col].abs()) / 2 + eps)
    ).mean()


def _compute_smdape(
    df: pd.DataFrame, model_col: str, rfc_col: str, eps: float
) -> float:
    """
    Calculate Symmetric Median Absolute Percentage Error (sMdAPE).

    sMdAPE provides skewness information to complement sMAPE.

    Formula: sMdAPE = Median[|ŷᵢ - yᵢ| / ((|ŷᵢ| + |yᵢ|)/2 + ε)]

    Reference
    ---------
    Wen et al. (2022), Applied Energy 325, 119906, Table 1
    """
    return (
        (df[model_col] - df[rfc_col]).abs()
        / ((df[model_col].abs() + df[rfc_col].abs()) / 2 + eps)
    ).median()


def _compute_rmsle(df: pd.DataFrame, model_col: str, rfc_col: str, eps: float) -> float:
    """
    Calculate Root Mean Square Logarithmic Error (RMSLE).

    RMSLE is complementary to percentage errors since it shows the logarithmic
    deviation values instead of percentage ones.

    Formula: RMSLE = √[(1/n) * Σ[log(1 + ((ŷᵢ + ε) - (yᵢ + ε))/(yᵢ + ε))]²]

    Reference
    ---------
    Wen et al. (2022), Applied Energy 325, 119906, Table 1
    """
    return np.sqrt(
        (
            np.log(
                1 + ((df[model_col] + eps) - (df[rfc_col] + eps)) / (df[rfc_col] + eps)
            )
            ** 2
        ).mean()
    )


def _compute_growth_error(
    df: pd.DataFrame, model_col: str, rfc_col: str, eps: float
) -> float:
    """
    Calculate Growth Error.

    The growth error shows the error in temporal scale. This indicator is ignored for dynamic time series.

    Formula: Growth error = ĝₜ - gₜ, where gₜ = [ln(yₜ) - ln(yₜ₀)] / (t - t₀)

    Reference
    ---------
    Wen et al. (2022), Applied Energy 325, 119906, Table 1
    """
    if len(df) < 2:
        logger.warning("Insufficient data for growth error calculation")
        return np.nan

    # Sort by time to ensure proper chronological order
    df_sorted = df.sort_values("year")

    # Validate time span
    t0 = df_sorted.index.get_level_values("year")[0]
    t1 = df_sorted.index.get_level_values("year")[-1]
    time_span = t1 - t0

    if time_span == 0:
        logger.warning("Zero time span for growth error calculation")
        return np.nan

    # Calculate growth rates
    def _compute_growth_rate(values: pd.Series) -> float:
        y0 = max(values.iloc[0], eps)
        y1 = max(values.iloc[-1], eps)
        return (np.log(y1) - np.log(y0)) / time_span

    model_growth = _compute_growth_rate(df_sorted[model_col])
    rfc_growth = _compute_growth_rate(df_sorted[rfc_col])
    return model_growth - rfc_growth


def _compute_missing(df_na: pd.DataFrame) -> int:
    """
    Calculate missing carriers count.
    """
    return len(df_na.index.get_level_values("carrier").unique())


def compute_all_indicators(
    df: pd.DataFrame,
    table: str,
    model_col: str = "Open-TYNDP",
    rfc_col: str = "TYNDP 2024",
    eps: float = 1e-6,
    carrier: str = None,
    df_na: pd.DataFrame = None,
) -> pd.DataFrame:
    """
    Compute all accuracy indicators for a given dataset.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing data for indicator calculations.
    table : str
        Benchmark metric to compute.
    model_col : str, default "Open-TYNDP"
        Column name for model/projected values (ŷᵢ).
    rfc_col : str, default "TYNDP 2024"
        Column name for reference/actual values (yᵢ).
    eps: float, default 1e-6
        Small value used when the denominator is zero.
    carrier: str, default None
        Name of the carrier for indicator calculation. If None, calculates overall table indicator.
    df_na : pd.DataFrame, default None
        DataFrame with missing values for missing carrier calculation.

    Returns
    -------
    pd.DataFrame
        DataFrame containing indicators, at carrier level if specified.
    """
    expected_cols = [model_col, rfc_col]
    if not df.columns.tolist() == expected_cols:
        logger.warning(f"Expected columns {expected_cols}, got {df.columns.tolist()}")
        return pd.DataFrame()

    indicators = {
        "sMPE": _compute_smpe(df, model_col, rfc_col, eps),
        "sMAPE": _compute_smape(df, model_col, rfc_col, eps),
        "sMdAPE": _compute_smdape(df, model_col, rfc_col, eps),
        "RMSLE": _compute_rmsle(df, model_col, rfc_col, eps),
    }

    if "snapshot" not in df.index.names:
        indicators["Growth Error"] = _compute_growth_error(df, model_col, rfc_col, eps)
    elif not carrier:
        indicators["Growth Error"] = "NA"

    if df_na is not None:
        indicators["Missing"] = _compute_missing(df_na)

    if carrier:
        indicators = {(table, carrier): indicators}
    else:
        indicators = {table: indicators}

    return pd.DataFrame.from_dict(indicators, orient="index")


def compute_indicators(
    df_raw: pd.DataFrame,
    table: str,
    options,
    carrier_col: str = "carrier",
    precision: int = 2,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate accuracy indicators following Wen et al. (2022) methodology to assess model performance
    against reference data. The function expects paired columns representing workflow estimates
    and TYNDP 2024 baseline values. The function computes both per-carrier and overall indicators.

    Hourly time series from the rfc_col will be aggregated to match the temporal resolution of model_col.

    Computes six key accuracy indicators:
    - Missing: Count of carrier dropped due to missing values
    - sMPE: Symmetric Mean Percentage Error (directional)
    - sMAPE: Symmetric Mean Absolute Percentage Error (magnitude)
    - sMdAPE: Symmetric Median Absolute Percentage Error (skewness)
    - RMSLE: Root Mean Square Logarithmic Error (logarithmic deviations)
    - Growth_Error: Growth Error (temporal analysis of transition rates)

    Reference: Wen, X., Jaxa-Rozen, M., Trutnevyte, E., 2022. Accuracy indicators for evaluating
    retrospective performance of energy system models. Applied Energy 325, 119906.
    https://doi.org/10.1016/j.apenergy.2022.119906

    Parameters
    ----------
    df_raw : pd.DataFrame
        DataFrame containing data for indicator calculations.
    table : str
        Benchmark metric to compute.
    options : dict
        Full benchmarking configuration.
    carrier_col : str, default "carrier"
        Column name for carrier/technology grouping.
    precision: int, default 2
        Number of decimal places to round to.

    Returns
    -------
    pd.DataFrame
       DataFrame with per carriers accuracy indicators.
    pd.Series
       Series containing overall accuracy indicators.
    """
    opt = options["tables"][table]

    # Aggregate time-varying data to the given snapshots
    if opt["table_type"] == "time_series":
        df_agg = match_temporal_resolution(df_raw)
    else:
        df_agg = df_raw

    mask = df_agg.isna().any(axis=1)
    df = df_agg[~mask]
    df_na = df_agg[mask]

    # Compute overall indicators of the table
    indicators = compute_all_indicators(df, table, df_na=df_na).round(precision)

    # Compute per-carrier indicators
    df_carrier = [
        compute_all_indicators(df_c, table, carrier=carrier)
        for carrier, df_c in df.groupby(level=carrier_col)
    ]
    missing_carriers = set(df_na.index.get_level_values("carrier"))
    df_carrier.extend(
        [pd.DataFrame(index=[(table, carrier) for carrier in missing_carriers])]
    )
    df_carrier = pd.concat(df_carrier).round(precision)

    return df_carrier, indicators


def compare_sources(
    table: str, benchmarks_raw: pd.DataFrame, options: dict
) -> tuple[pd.DataFrame, pd.Series]:
    """
    Compare data sources for a specified table using accuracy indicators. The function expects
    paired columns representing workflow results estimates and TYNDP 2024 baseline values.

    Parameters
    ----------
    table : str
        Benchmark metric to compute.
    benchmarks_raw: pd.DataFrame
        Combined DataFrame containing both Open-TYNDP and TYNDP 2024 data.
    options : dict
        Full benchmarking configuration.

    Returns
    -------
    pd.DataFrame
       DataFrame containing original data with appended multi-value accuracy metric columns.
    pd.Series
       Series containing single-value accuracy metrics.
    """

    # Clean data
    logger.info(f"Making benchmark for {table} using TYNDP 2024 and Open-TYNDP")
    benchmarks = benchmarks_raw.query("table==@table").dropna(how="all", axis=1)
    available_columns = [
        c for c in benchmarks.columns if c not in ["value", "source", "unit"]
    ]
    df = benchmarks.pivot_table(
        index=available_columns, values="value", columns="source", dropna=False
    )

    # Compare sources
    df, indicators = compute_indicators(df, table, options)

    return df, indicators


def compute_overall_accuracy(
    benchmarks_raw: pd.DataFrame, options: dict
) -> pd.DataFrame:
    """
    Compute overall accuracy indicators for a specified benchmark table.

    Parameters
    ----------
    benchmarks_raw: pd.DataFrame
        Combined DataFrame containing both Open-TYNDP and TYNDP 2024 data.
    options : dict
        Full benchmarking configuration.

    Returns
    -------
    pd.Series
       Series containing overall accuracy metrics.
    """
    logger.info("Making global benchmark using TYNDP 2024 and Open-TYNDP")
    tables_series = [  # noqa: F841
        t for t, v in options["tables"].items() if v["table_type"] == "time_series"
    ]
    df_global = (
        benchmarks_raw.query("table not in @tables_series")
        .dropna(how="all", axis=1)
        .pivot_table(
            index=["scenario", "year", "carrier"], values="value", columns="source"
        )
    )
    mask = df_global.isna().any(axis=1)
    df = df_global[~mask]
    df_na = df_global[mask]
    indicator_total = compute_all_indicators(
        df, "Total (excl. time series)", df_na=df_na
    ).round(2)

    return indicator_total


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "make_benchmark",
            opts="",
            clusters="all",
            sector_opts="",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    options = snakemake.params["benchmarking"]
    scenario = "TYNDP " + snakemake.params["scenario"]
    benchmarks_fn = snakemake.input.benchmarks
    results_fn = snakemake.input.results

    # Load data
    benchmarks_raw = load_data(benchmarks_fn, results_fn, scenario)

    # Compute benchmarks
    logger.info("Computing benchmarks")

    tqdm_kwargs = {
        "ascii": False,
        "unit": " benchmark",
        "total": len(options["tables"]),
        "desc": "Computing benchmark",
    }

    func = partial(
        compare_sources,
        benchmarks_raw=benchmarks_raw,
        options=options,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        results = list(tqdm(pool.imap(func, options["tables"].keys()), **tqdm_kwargs))
        benchmarks, indicators = zip(*results)

    # Combine and write all benchmark data
    os.makedirs(snakemake.output.benchmarks, exist_ok=True)
    for benchmark in benchmarks:
        if not benchmark.empty:
            table = benchmark.index.get_level_values(0)[0]
            benchmark.loc[table].to_csv(
                snakemake.output.benchmarks
                + f"/{table}_s_{snakemake.wildcards.clusters}_{snakemake.wildcards.opts}_{snakemake.wildcards.sector_opts}_all_years.csv"
            )

    # Compute global indicator
    indicators_total = compute_overall_accuracy(benchmarks_raw, options)

    indicators = pd.concat(list(indicators) + [indicators_total])
    indicators.to_csv(snakemake.output.kpis)
