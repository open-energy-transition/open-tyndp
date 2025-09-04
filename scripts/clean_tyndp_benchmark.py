# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script cleans the raw TYNDP 2024 Scenarios Report Data that will be used for benchmarking.

It reads and processes all the tables defined in the configuration file. The correct indexes
and headers are then assigned. The data structure is subsequently converted to a long format.
Finally, the units are converted to standard units: MW for power units and MWh for energy units.
"""

import logging
import multiprocessing as mp
from functools import partial

import pandas as pd
from tqdm import tqdm

from scripts._helpers import SCENARIO_DICT, configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def _safe_sheet(sn, scenario):
    if isinstance(sn, dict):
        return sn[scenario]
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
    for hdr in header:
        header_array.append(df.iloc[hdr - 1].infer_objects(copy=False).ffill().values)

    if not header_array:
        return df

    if len(header_array) == 1:
        df.columns = header_array[0]
    else:
        df.columns = pd.MultiIndex.from_arrays(header_array)

    df = df.iloc[max(header) :]
    df = df.dropna(how="all").infer_objects(copy=False).fillna(0.0)

    df.columns.names = names

    return df


def _convert_units(
    df: pd.DataFrame,
    source_unit: str,
    unit_conversion: dict[str, float],
    value_col: str = "value",
) -> pd.DataFrame:
    """
    Convert units and add unit columns.

    Automatically determines target unit based on source unit type:
    - Energy units (TWh, GWh, MWh, kWh) → MWh
    - Power units (GW, MW, kW) → MW

    Parameters
    ----------
    df : pd.DataFrame
        Long-format DataFrame containing values to convert.
    source_unit : str
        Source unit of the values.
    unit_conversion : dict[str, float]
        Dictionary mapping units to conversion factors (to base unit).
    value_col : str, default "value"
        Name of the column containing values to convert.

    Returns
    -------
    pd.DataFrame
        DataFrame with converted values and unit columns added.
    """
    # Determine target unit based on source unit type
    energy_units = {"TWh", "GWh", "MWh", "kWh"}
    power_units = {"GW", "MW", "kW"}

    if source_unit in energy_units:
        target_unit = "MWh"
    elif source_unit in power_units:
        target_unit = "MW"
    else:
        # Unknown unit type, keep original
        target_unit = source_unit

    # Convert values using conversion factor from config
    conversion_factor = unit_conversion[source_unit]
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce") * conversion_factor
    df["unit"] = target_unit

    return df


def _add_identifier(s):
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
    if table_type not in ["scenario_comparison", "time_serie"]:
        logger.warning(f"Table type '{table_type}' not implemented yet")
        return pd.DataFrame()

    # Parameters
    sheet_name = _safe_sheet(opt["sheet_name"], scenario)
    df = benchmarks_raw[sheet_name]
    table_config = options["table_types"][opt["table_type"]]
    nrows = opt.get("nrows", None)
    ncolumns = opt.get("ncolumns", None)
    names = opt["names"]

    index_col = opt.get("index_col", table_config["index_col"])
    if isinstance(index_col, int):
        index_col = [index_col]

    header = opt.get("header", table_config["header"])
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
    source_unit = opt["unit"]
    unit_conversion = options["processing"]["unit_conversion"]
    df_converted = _convert_units(df_long, source_unit, unit_conversion)

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
    df_converted = df_converted[
        ~df_converted.carrier.isin(["sum", "aggregated", "total generation", "total"])
    ]

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
        _safe_sheet(j["sheet_name"], scenario) for i, j in options["tables"].items()
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
