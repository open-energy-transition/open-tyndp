# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script cleans the raw TYNDP 2024 Scenarios Report Data that will be used for benchmarking.

It reads and processes all the tables defined in the configuration file. The correct indexes
and headers are then assigned. The data structure is subsequently converted to a long format.
Finally, the units are converted to standard units: MW for power units and MWh for energy units.
"""

import datetime
import logging
import multiprocessing as mp
from functools import partial

import pandas as pd
from tqdm import tqdm

from scripts._helpers import (
    SCENARIO_DICT,
    configure_logging,
    convert_units,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def _safe_sheet(sn: str | dict, scenario: str) -> str:
    if isinstance(sn, dict):
        return sn.get(scenario, "")
    return sn


def _process_index(
    df: pd.DataFrame, index_col: list[int], names: list[str], nrows: int = None
) -> pd.DataFrame:
    """
    Process index columns to create proper index structure.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with indices in data columns.
    index_col : list[int]
        List of Excel column indices for index.
    names : list[str]
        List of index names.
    nrows : int, optional
        Number of rows to process. If None (default), process all rows.

    Returns
    -------
    pd.DataFrame
        DataFrame with proper index.
    """
    if nrows:
        df = df[:nrows]

    index_array = [i - 1 for i in index_col]
    df = df.iloc[:, min(index_array) :]
    if len(index_col) > 1:
        df.loc[:, index_array] = df.loc[:, index_array].ffill()

    df = df.set_index(index_array)
    df.index.names = names

    return df


def _process_header(
    df: pd.DataFrame, header: list[int], names: list[str], ncolumns: int = None
) -> pd.DataFrame:
    """
    Process header rows to create proper column structure.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with headers in data rows.
    header : list[int]
        List of Excel row indices for headers.
    names : list[str]
        List of column names.
    ncolumns : int, optional
        Number of columns to process. If None (default), use all columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with proper column headers.
    """
    if ncolumns:
        df = df.iloc[:, : (ncolumns - df.index.nlevels) - 1]

    header_array = []
    for hdr_i in header:
        hdr = df.iloc[hdr_i - 1]
        if len(header) > 1:
            hdr = hdr.infer_objects(copy=False).ffill()
        header_array.append(hdr.values)

    if not header_array:
        return df

    if len(header_array) == 1:
        df.columns = header_array[0]
    else:
        df.columns = pd.MultiIndex.from_arrays(header_array)

    df = df.iloc[max(header) :].dropna(how="all")
    df = df.loc[df.index.dropna(), df.columns.dropna()]
    df = df.infer_objects(copy=False).fillna(0.0)

    df.columns.names = names

    return df


def _add_identifier(s: str) -> str:
    """
    Add institution identifier to scenario name.
    """
    if s.startswith("TYNDP") or s.startswith("EC"):
        return s
    elif any(i in s for i in ["DE", "GA", "NT"]):
        return "TYNDP " + s
    elif any(i in s for i in ["IA", "S3"]):
        return "EC IA S3"
    else:
        return s


def load_benchmark(
    benchmarks_raw: dict[str, pd.DataFrame], table: str, scenario: str, options: dict
) -> pd.DataFrame:
    """
    Load and process benchmark data from TYNDP Excel sheets.

    Parameters
    ----------
    benchmarks_raw : dict[str, pd.DataFrame]
        Dictionary of pandas DataFrames from Excel sheets, keyed by sheet name.
    table : str
        Name of table to load.
    scenario: str
        Name of scenario to load.
    options : dict
        Full benchmarking configuration.

    Returns
    -------
    pd.DataFrame
        Cleaned benchmark data in long format.
    """

    # Handle each table type
    opt = options["tables"][table]
    table_type = opt["table_type"]
    if table_type not in ["scenario_comparison", "time_series"]:
        logger.warning(f"Table type '{table_type}' not implemented yet")
        return pd.DataFrame()

    # Parameters
    sheet_name = _safe_sheet(opt["sheet_name"], scenario)
    if sheet_name == "":
        logger.warning(f"No sheet name found for {table} in {scenario}")
        return pd.DataFrame()
    df = benchmarks_raw[sheet_name]
    nrows = opt.get("nrows", None)
    ncolumns = opt.get("ncolumns", None)
    names = opt["names"]

    # Fix temporal labeling - source data uses 00:00 as end-of-period (previous day's last hour);
    # convert to beginning-of-period (current day's first hour)
    if table == "generation_profiles":
        time_col = df.iloc[
            5:, 0
        ]  # Start reading from row 6 where actual snapshot data begins
        datetime_series = pd.to_datetime(time_col)
        midnight_mask = datetime_series.dt.time == datetime.time(0, 0)
        df.loc[time_col.index[midnight_mask], df.columns[0]] = datetime_series[
            midnight_mask
        ] + pd.Timedelta(days=1)

    index_col = opt["index_col"]
    if isinstance(index_col, int):
        index_col = [index_col]

    header = opt["header"]
    if isinstance(header, int):
        header = [header]

    # Process row indexes and headers
    df = _process_index(df, index_col, names["index"], nrows)
    df = _process_header(df, header, names["column"], ncolumns)

    # Convert to long format
    if df.columns.nlevels > df.index.nlevels:
        id_vars = [
            tuple(df.index.names + [""] * (df.columns.nlevels - df.index.nlevels))
        ]
    else:
        id_vars = df.index.names
    df_long = df.reset_index().melt(id_vars=id_vars)
    df_long.columns = [i[0] if isinstance(i, tuple) else i for i in df_long.columns]

    # Apply unit conversion
    df_long["unit"] = opt["unit"]
    df_converted = convert_units(df_long)

    # Add table identifier
    df_converted["table"] = table

    # Apply exceptions
    if table == "generation_profiles":
        # TODO Validate planning year assumption
        df_converted["year"] = 2040
        df_converted["scenario"] = scenario

    # Clean data
    if "scenario" in df_converted.columns:
        df_converted["scenario"] = (
            df_converted["scenario"]
            .replace(SCENARIO_DICT, regex=True)
            .apply(_add_identifier)
        )

    df_converted["year"] = df_converted["year"].astype(int)
    df_converted["carrier"] = df_converted["carrier"].str.lower().str.rstrip("* ")

    # Remove aggregated values except for electricity demand
    if table != "elec_demand":
        df_converted = df_converted[
            ~df_converted.carrier.isin(
                ["sum", "aggregated", "total generation", "total"]
            )
        ]
    # Keep aggregated electricity demand - Input data for Open-TYNDP is provided in aggregated form
    else:
        df_converted = df_converted[df_converted.carrier == "aggregated"]

    return df_converted


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("clean_tyndp_benchmark")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    options = snakemake.params["benchmarking"]
    scenario = snakemake.params["scenario"]

    # Read benchmarks
    logger.info("Reading raw benchmark data")
    sheet_names = [
        _safe_sheet(j["sheet_name"], scenario)
        for i, j in options["tables"].items()
        if _safe_sheet(j["sheet_name"], scenario)
    ]
    benchmarks_raw = pd.read_excel(
        snakemake.input.scenarios_figures, sheet_name=sheet_names, header=None
    )

    logger.info("Parsing benchmark data")
    tqdm_kwargs = {
        "ascii": False,
        "unit": " benchmark",
        "total": len(options["tables"]),
        "desc": "Loading TYNDP benchmark data",
    }

    func = partial(
        load_benchmark,
        benchmarks_raw,
        scenario=scenario,
        options=options,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        benchmarks = list(
            tqdm(pool.imap(func, options["tables"].keys()), **tqdm_kwargs)
        )

    # Combine all benchmark data
    benchmarks_combined = pd.concat(benchmarks, ignore_index=True).assign(
        source="TYNDP 2024"
    )
    if benchmarks_combined.empty:
        logger.warning("No benchmark data was successfully processed")

    # Save data
    benchmarks_combined.to_csv(snakemake.output.benchmarks, index=False)
