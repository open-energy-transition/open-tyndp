# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
import re
import socket
from bisect import bisect_right
from collections.abc import Callable

import numpy as np
import pandas as pd

from scripts._helpers import logger, make_index

SCENARIO_DICT = {
    "Distributed Energy": "DE",
    "Global Ambition": "GA",
    r"National Trends\s*\+": "NT",
    r"NT\s*\+": "NT",
    "National Trends": "NT",
}
ENERGY_UNITS = {"TWh", "GWh", "MWh", "kWh"}
POWER_UNITS = {"GW", "MW", "kW"}
PRICE_UNITS = {"EUR/MWh", "EUR/MWh_e", "EUR/MWh_H2"}


def extract_grid_data_tyndp(
    links,
    replace_dict: dict = {},
    expand_from_index: bool = True,
    idx_prefix: str = "",
    idx_connector: str = "",
    idx_suffix: str = "",
    idx_separator: str = " ",
):
    """
    Extract TYNDP reference grid data from the raw input table.

    Parameters
    ----------
    links : pd.DataFrame
        DataFrame with raw links to extract grid information from.
    replace_dict : dict
        Dictionary with region names to replace.
    expand_from_index : bool
        Whether to expand the bus0 and bus1 from index or directly use the columns.
    idx_prefix : str, optional
        Prefix to prepend to generated indices.
    idx_connector : str, optional
        Separator string between bus0 and bus1 in generated indices (e.g., "->", "-").
    idx_suffix : str, optional
        Suffix to append to generated indices.
    idx_separator : str, optional
        String used to join index components (prefix, bus0, connector, bus1, suffix).

    Returns
    -------
    pd.DataFrame
        DataFrame with extracted grid data information with nominal capacity in input unit, bus0 and bus1.
    """

    if expand_from_index:
        links.loc[:, "Border"] = links["Border"].replace(replace_dict, regex=True)
        links = pd.concat(
            [
                links,
                links.Border.str.split("-", expand=True).set_axis(
                    ["bus0", "bus1"], axis=1
                ),
            ],
            axis=1,
        )
    elif "bus0" in links.columns and "bus1" in links.columns:
        links.loc[:, ["bus0", "bus1"]] = links[["bus0", "bus1"]].replace(
            replace_dict, regex=True
        )
    else:
        raise KeyError(
            f"Columns 'bus0' and 'bus1' must be present in the input DataFrame. Available columns: {list(links.columns)}"
        )

    # Create forward and reverse direction dataframes
    forward_links = links[["bus0", "bus1", "Summary Direction 1"]].rename(
        columns={"Summary Direction 1": "p_nom"}
    )

    reverse_links = links[["bus1", "bus0", "Summary Direction 2"]].rename(
        columns={"bus1": "bus0", "bus0": "bus1", "Summary Direction 2": "p_nom"}
    )

    # Combine into unidirectional links and return
    links = pd.concat([forward_links, reverse_links])

    links.index = links.apply(
        make_index,
        axis=1,
        prefix=idx_prefix,
        connector=idx_connector,
        suffix=idx_suffix,
        separator=idx_separator,
    )

    return links


def safe_pyear(
    year: int | str,
    available_years: list[int] = [2030, 2040, 2050],
    source: str = "TYNDP",
    verbose: bool = True,
) -> int:
    """
    Checks and adjusts whether a given pyear is in the available years of a given data source. If not, it
    falls back to the previous available year.

    Parameters
    ----------
    year : int
        Planning horizon year which will be checked and possibly adjusted to previous available year.
    available_years : list[int], optional
        List of available years. Default is [2030, 2040, 2050].
    source : str, optional
        Source of the data for which availability will be checked. For logging purpose only. Default is "TYNDP".
    verbose : bool, optional
        Whether to activate verbose logging. Default is True.

    Returns
    -------
    year_new : int
        Safe pyear adjusted for available years.
    """

    if not available_years:
        raise ValueError(
            "No `available_years` provided. Expected a non-empty list of years."
        )
    if not isinstance(year, int):
        year = int(year)
    if year not in available_years:
        year_new = available_years[
            bisect_right(sorted(available_years), year, lo=1) - 1
        ]
        if verbose:
            logger.warning(
                f"{source} data unavailable for planning horizon {year}. Falling back to previous available year {year_new}."
            )
    else:
        year_new = year

    return year_new


def map_tyndp_carrier_names(
    df: pd.DataFrame,
    carrier_mapping_fn: str,
    on_columns: list[str],
    drop_on_columns=False,
):
    """
    Map external carriers to available tyndp_carrier names based on an input mapping. Optionally drop merged on columns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with external carriers to map
    carrier_mapping_fn : str
        Path to file with mapping from external carriers to available tyndp_carrier names.
    on_columns : list[str]
        Columns to merge on between the external carriers and tyndp_carriers.
    drop_on_columns : bool, optional
        Whether to drop merge columns and rename `open_tyndp_carrier` and `open_tyndp_index` to `carrier`
        and `index_carrier`. Default is False.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with external carriers mapped to available tyndp_carriers and index_carriers.
    """

    # Read TYNDP carrier mapping
    carrier_mapping = (
        pd.read_csv(carrier_mapping_fn)[
            on_columns + ["open_tyndp_carrier", "open_tyndp_index", "open_tyndp_type"]
        ]
    ).dropna()

    # Map the carriers
    df = df.merge(carrier_mapping, on=on_columns, how="left")

    # If the carrier is DSR or Other Non-RES, the different price bands are too diverse for a robust external
    # mapping. Instead, we will combine the carrier and type information.
    if "pemmdb_carrier" in on_columns:

        def normalize_carrier(s):
            return s.lower().replace(" ", "-").replace("other-non-res", "chp")

        # Other Non-RES are assumed to represent CHP plants (according to TYNDP 2024 Methodology report p.37)
        df = df.assign(
            open_tyndp_carrier=lambda x: np.where(
                x["pemmdb_carrier"].isin(["DSR", "Other Non-RES"]),
                x["pemmdb_carrier"].apply(normalize_carrier),
                x["open_tyndp_carrier"],
            ),
            open_tyndp_index=lambda x: np.where(
                x["pemmdb_carrier"].isin(["DSR", "Other Non-RES"]),
                x["open_tyndp_carrier"]
                + "-"
                + x["pemmdb_type"].apply(normalize_carrier),
                x["open_tyndp_index"],
            ),
        )

    if not drop_on_columns:
        return df

    # Otherwise drop merge columns and rename to new "carrier" and "index_carrier" column
    df = df.drop(on_columns, axis="columns").rename(
        columns={
            "open_tyndp_carrier": "carrier",
            "open_tyndp_index": "index_carrier",
        }
    )

    # Move "carrier" and "index_carrier" to the front
    cols = ["carrier", "index_carrier"] + [
        col for col in df.columns if col not in ["carrier", "index_carrier"]
    ]

    return df[cols]


def convert_units(
    df: pd.DataFrame,
    unit_col: str = "unit",
    value_col: str = "value",
    invert: bool = False,
) -> pd.DataFrame:
    """
    Convert values to standardized units based on unit type.

    Energy values are converted to MWh, power values to MW:
    - Energy units (TWh, GWh, MWh, kWh) → MWh
    - Power units (GW, MW, kW) → MW

    When invert=False (default):
        - Values are converted from unit_col units to standard units (MWh/MW)
        - The "unit" column is updated to reflect the standardized unit

    When invert=True:
        - Values are converted from standard units (MWh/MW) back to unit_col units
        - The "unit" column is NOT modified
        - Useful for reverting previously standardized data

    Parameters
    ----------
    df : pd.DataFrame
        Long-format DataFrame containing values to convert.
    unit_col : str, default "unit"
        Name of the column containing the unit information.
        When invert=False: contains source units to convert from.
        When invert=True: contains target units to convert to.
    value_col : str, default "value"
        Name of the column containing values to convert.
    invert : bool, default False
        If False, convert to standard units and update "unit" column.
        If True, convert from standard units using inverse factors without modifying "unit" column.

    Returns
    -------
    pd.DataFrame
        DataFrame with converted values.
    """
    df = df.copy()

    unit_conversion = {
        "TWh": 1000000,
        "GWh": 1000,
        "MWh": 1,
        "GW": 1000,
        "MW": 1,
        "kW": 0.001,
    }

    if invert:
        # Inverse conversion factor to revert unit
        unit_conversion = {k: 1 / v for k, v in unit_conversion.items()}

    # Convert values using conversion factors
    conversion_factors = df[unit_col].map(unit_conversion)
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce") * conversion_factors

    # Update unit column
    if not invert:
        df["unit"] = df[unit_col].apply(
            lambda x: "MWh" if x in ENERGY_UNITS else "MW" if x in POWER_UNITS else x
        )

    return df


def check_cyear(cyear: int, scenario: str) -> int:
    """
    Check if the climatic year is valid for the given scenario.

    Parameters
    ----------
    cyear : int
        Climatic year to validate.
    scenario : str
        TYNDP scenario name.

    Returns
    -------
    int
        Valid climatic year, falling back to 2009 if the input is not available.
    """

    valid_years = {
        "NT": [1995, 2008, 2009],
        "DE": [1995, 2008, 2009],
        "GA": [1995, 2008, 2009],
    }

    if cyear not in valid_years[scenario]:
        logger.warning(
            f"Snapshot year {cyear} doesn't match available TYNDP data. Falling back to 2009."
        )
        cyear = 2009

    return cyear


def get_tyndp_conventional_thermals(
    mapping: pd.DataFrame,
    tyndp_conventional_carriers: list[str],
    group_conventionals: bool,
    include_h2_fuel_cell: bool,
    include_h2_turbine: bool,
) -> tuple[dict[str, str], list[str]]:
    """
    Get list of TYNDP conventional thermal generation technologies.

    Parameters
    ----------
    mapping : pd.DataFrame
        TYNDP conventional carrier mapping (grouped or ungrouped).
    tyndp_conventional_carriers : list[str]
        TYNDP conventional carriers.
    group_conventionals : bool
        Whether to group conventional thermal technologies.
    include_h2_fuel_cell : bool
        Whether to include hydrogen fuel cell technology.
    include_h2_turbine : bool
        Whether to include hydrogen turbine technology.

    Returns
    -------
    tuple[dict[str, str], list[str]]
        Dictionary with conventional mapping and list of conventional thermal technology names.
    """

    # Filter mapping for conventional carriers while setting oil as common fuel for oil technologies
    mapping = (
        mapping[["open_tyndp_carrier", "open_tyndp_type", "pypsa_eur_carrier"]]
        .query("open_tyndp_carrier in @tyndp_conventional_carriers")
        .loc[lambda df: df.index.notna()]
        .replace(
            {"open_tyndp_carrier": ["oil-light", "oil-heavy", "oil-shale"]}, "oil"
        )  # TODO To remove once the three carriers have been implemented
    )

    if group_conventionals:
        mapping = mapping.groupby("open_tyndp_type").first()

    conventional_dict = mapping.open_tyndp_carrier.to_dict()
    conventional_thermals = list(conventional_dict)

    if include_h2_fuel_cell:
        conventional_thermals.append("h2-fuel-cell")
    if include_h2_turbine:
        conventional_thermals.append("h2-ccgt")

    return conventional_dict, conventional_thermals


def interpolate_demand(
    available_years: list[int],
    pyear: int,
    load_single_year_func: Callable,
    **load_kwargs,
) -> pd.DataFrame | pd.Series:
    """
    Interpolate demand between available years.

    Parameters
    ----------
    available_years : list[int]
        Sorted list of years for which data is available.
    pyear : int
        Planning year to interpolate demand for.
    load_single_year_func : Callable
        Function to load data for a single planning year.
    **load_kwargs
        Keyword arguments to pass to load_single_year_func. Must include 'pyear'
        as a parameter key, which will be overridden with interpolation boundary years.

    Returns
    -------
    pd.DataFrame | pd.Series
        Interpolated demand data.
    """
    # Currently, only interpolation is implemented, not extrapolation
    idx = bisect_right(available_years, pyear)
    if idx == 0:
        # Planning horizon is before all available years
        logger.warning(
            f"Year {pyear} is before the first available year {available_years[0]}. "
            f"Falling back to first available year."
        )
        year_lower = year_upper = available_years[0]
    elif idx == len(available_years):
        # Planning horizon is after all available years
        logger.warning(
            f"Year {pyear} is after the latest available year {available_years[-1]}. "
            f"Falling back to latest available year."
        )
        year_lower = year_upper = available_years[-1]
    else:
        year_lower = available_years[idx - 1]
        year_upper = available_years[idx]

    logger.debug(f"Interpolating {pyear} from {year_lower} and {year_upper}")

    kwargs_lower = {**load_kwargs, "pyear": year_lower}
    kwargs_upper = {**load_kwargs, "pyear": year_upper}

    df_lower = load_single_year_func(**kwargs_lower)
    df_upper = load_single_year_func(**kwargs_upper)

    # Check if data was loaded successfully
    if df_lower.empty and df_upper.empty:
        logger.error("Both years failed to load")
        return pd.DataFrame()
    elif df_lower.empty:
        logger.warning(
            f"Year {year_lower} failed to load. Filling with zeros for interpolation."
        )
        df_lower = pd.DataFrame(0, index=df_upper.index, columns=df_upper.columns)
    elif df_upper.empty:
        logger.warning(
            f"Year {year_upper} failed to load. Using data from lower year for interpolation."
        )
        df_upper = df_lower

    if year_upper == year_lower:
        return df_lower

    # Handle column mismatches for DataFrames (only relevant for DataFrame, not Series)
    if isinstance(df_lower, pd.DataFrame) and isinstance(df_upper, pd.DataFrame):
        missing_in_lower = df_upper.columns.difference(df_lower.columns)
        missing_in_upper = df_lower.columns.difference(df_upper.columns)

        if len(missing_in_lower) > 0 or len(missing_in_upper) > 0:
            logger.warning(
                f"Column mismatch between {year_lower} and {year_upper}. "
                f"Missing columns filled with zeros. "
                f"Missing in {year_lower}: {list(missing_in_lower)}, "
                f"Missing in {year_upper}: {list(missing_in_upper)}"
            )
        df_lower_aligned, df_upper_aligned = df_lower.align(
            df_upper, join="outer", axis=1, fill_value=0
        )
    else:
        # For Series, just align
        df_lower_aligned, df_upper_aligned = df_lower.align(
            df_upper, join="outer", fill_value=0
        )

    # Perform linear interpolation
    weight = (pyear - year_lower) / (year_upper - year_lower)
    result = df_lower_aligned * (1 - weight) + df_upper_aligned * weight

    return result


def find_free_port(start_port: int = 8050, max_attempts: int = 50) -> int:
    """
    Find the first available port starting from start_port.

    Parameters
    ----------
    start_port : int, optional
        Port number to begin scanning from. Default is 8050.
    max_attempts : int, optional
        Maximum number of ports to check before raising an error. Default is 50.

    Returns
    -------
    int
        First available port number in the scanned range.
    """
    for port in range(start_port, start_port + max_attempts):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(("127.0.0.1", port))
                return port
        except OSError:
            continue
    raise RuntimeError(
        f"Could not find free port in range {start_port}-{start_port + max_attempts}"
    )


def _add_new_profiles_to_existing(
    component_t: dict,
    attr: str,
    new_profiles: pd.DataFrame,
) -> None:
    """
    Safely add new time series data to a PyPSA network component.

    Parameters
    ----------
    component_t : dict
        The time-varying PyPSA component (e.g., n.generators_t).
    attr : str
        The attribute name (e.g., 'p_max_pu', 'inflow').
    new_profiles : pd.DataFrame
        New data to concatenate.
    """
    if new_profiles.empty:
        return

    existing = component_t[attr]
    index_name = existing.index.name

    component_t[attr] = pd.concat([existing, new_profiles], axis=1)
    component_t[attr].index.name = index_name


def align_demand_to_snapshots(
    demand: pd.DataFrame, snapshots: pd.DatetimeIndex, format: str = None
) -> pd.DataFrame:
    """
    Convert demand index to DatetimeIndex, adjust year to match snapshots, and reindex to snapshots.

    Parameters
    ----------
    demand : pd.DataFrame
        Demand time series with a datetime-compatible index.
    snapshots : pd.DatetimeIndex
        Target snapshot index to align demand to.
    format : str, optional
        Datetime format string for parsing the demand index. Default is None.

    Returns
    -------
    pd.DataFrame
        Demand data reindexed to the provided snapshots.
    """

    demand.index = pd.to_datetime(demand.index, format=format)
    target_year = snapshots[0].year
    demand.index = demand.index.map(lambda x: x.replace(year=target_year))

    return demand.reindex(snapshots)


def extract_crossborder_pattern(df: pd.DataFrame, connector: str = "-"):
    return df.index.str.extract(
        rf"^(.*?)(\w+)( H2)*{re.escape(connector)}(\w+)( H2)*(.*)$"
    ).fillna("")


def normalize_direction(
    df: pd.DataFrame | pd.Series,
    cols: list[str] | None = None,
    buses_from_index: bool = False,
    connector: str = "->",
    format_index: bool = False,
) -> pd.DataFrame:
    """
    Normalize link direction so bus0 is always alphabetically before bus1.

    Parameters
    ----------
    df : pd.DataFrame | pd.Series
        Input data to normalize.
    cols : list[str] | None, default=None
        Names of flow columns whose sign should be flipped when direction is
        swapped. Required when `df` is a DataFrame and ignored when `df` is a
        Series.
    buses_from_index : bool, default=False
        If True, parse index instead of reading
        "bus0"/"bus1" columns directly.
    connector : str, default="->"
        String separator between bus0 and bus1 in the index.
    format_index : bool, default=False
        If True, reformat the index after normalization to
        ``"bus0[_suffix]->bus1[_suffix]"``.

    Returns
    -------
    pd.DataFrame | pd.Series
        Data with normalized link direction.
    """
    is_series = isinstance(df, pd.Series)

    if is_series:
        df = df.to_frame("value")
        cols = ["value"]
    elif isinstance(df, pd.DataFrame):
        if cols is None or len(cols) == 0:
            raise ValueError(
                "No column names `cols` to normalize. `cols` must be provided when `df` is a DataFrame."
            )
    else:
        raise TypeError("`df` must be a pandas Series or DataFrame.")

    if buses_from_index:
        df.loc[
            :, ["prefix", "bus0", "bus0_suffix", "bus1", "bus1_suffix", "suffix"]
        ] = extract_crossborder_pattern(df, connector).values
    # Protect swapping of direction for H2 import links and the bus0 starting with "X"
    mask_import_link = df["prefix"].str.contains("H2 import", na=False)
    mask = (
        (df["bus0"] > df["bus1"]) & ~df["bus0"].str.startswith("X") & ~mask_import_link
    )
    assignments = {col: np.where(mask, -df[col], df[col]) for col in cols}

    # Add bus0, bus1, and border to assignments
    assignments.update(
        {
            "bus0": np.where(mask, df["bus1"], df["bus0"]),
            "bus1": np.where(mask, df["bus0"], df["bus1"]),
            "border": lambda df: (
                df.prefix
                + df.bus0
                + df.bus0_suffix
                + connector
                + df.bus1
                + df.bus1_suffix
                + df.suffix
            ),
        }
    )

    df = df.assign(**assignments).drop(
        columns=["prefix", "bus0_suffix", "bus1_suffix", "suffix"], errors="ignore"
    )

    if buses_from_index:
        df = df.set_index("border")

    if format_index:
        idx_groups = extract_crossborder_pattern(df, connector)

        mask_h2_pipeline = idx_groups[0] == "H2 pipeline "
        idx_groups.loc[mask_h2_pipeline, [2, 4]] = " H2"

        mask_import = idx_groups[0].str.contains("H2 import")
        idx_groups.loc[mask_import, 1] = "X" + idx_groups.loc[mask_import, 1]
        idx_groups.loc[mask_import, 4] = " H2"

        df.index = idx_groups[1] + idx_groups[2] + "->" + idx_groups[3] + idx_groups[4]

    if buses_from_index:
        df.index.name = "border"

    if is_series:
        df = df.value

    return df
