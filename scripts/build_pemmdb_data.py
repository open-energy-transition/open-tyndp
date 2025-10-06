# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Collects and bundles the available PEMMDB v2.4 capacities and profiles for different PEMMDB technologies from TYNDP data bundle for a given planning horizon.

Outputs
-------
Cleaned CSV file with all NT capacities (p_nom) in long format and NetCDF file containing the must run obligations (p_min_pu) and availability (p_max_pu) for each of the different PEMMDB technologies.

- ``resources/pemmdb_capacities_{planning_horizon}.csv`` in long format
- ``resources/pemmdb_profiles_{planning_horizon}.nc`` with the following structure:

    ===================  ====================  =========================================================
    Field                Dimensions            Description
    ===================  ====================  =========================================================
    profile              time, bus, carrier,   the per unit hourly availability and must-run obligations
                         type                  for each bus and PEMMDB technology
    ===================  ====================  =========================================================
"""

import logging
import multiprocessing as mp
import sys
from functools import partial
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from scripts._helpers import (
    configure_logging,
    convert_units,
    get_snapshots,
    map_tyndp_carrier_names,
    safe_pyear,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

CONVENTIONALS = [
    "Nuclear",
    "Hard coal",
    "Lignite",
    "Gas",
    "Light oil",
    "Heavy oil",
    "Oil shale",
]

RENEWABLES = [
    "Solar",
    "Wind",
    "Hydro",
]

PEMMDB_SHEET_MAPPING = {
    "Gas": "Thermal",
    "Nuclear": "Thermal",
    "Hard coal": "Thermal",
    "Lignite": "Thermal",
    "Light oil": "Thermal",
    "Heavy oil": "Thermal",
    "Oil shale": "Thermal",
    "Hydrogen": "Thermal",
}


def read_pemmdb_data(
    node: str,
    pemmdb_dir: str,
    cyear: int,
    pyear: int,
    required_techs: list[str] = None,
) -> dict[str, dict[str, pd.DataFrame]]:
    """
    Read raw data from the PEMMDB for a specific planning and climate year, and a given set of technologies.

    Parameters
    ----------
    node : str
        Node name to read data for.
    pemmdb_dir : str
        Path to directory containing PEMMDB data.
    cyear : int
        Climate year to read data for.
    pyear : int
        Planning year used for data retrieval (fallback year if pyear_i not available).
    required_techs : list[str], optional
        List of required technologies to read PEMMDB data for.

    Returns
    -------
    dict[str, dict[str, pd.DataFrame]]
        A dictionary containing the read data for all given nodes in the form of {node: data}.
    """
    fn = Path(
        pemmdb_dir,
        str(pyear),
        f"PEMMDB_{node.replace('GB', 'UK')}_NationalTrends_{pyear}.xlsx",
    )

    if not fn.is_file():
        logger.info(f"No PEMMDB data available for {node} in {pyear}.")
        return None

    try:
        if required_techs:
            required_sheets = [
                PEMMDB_SHEET_MAPPING.get(tech, tech)
                for tech in required_techs
                if PEMMDB_SHEET_MAPPING.get(tech, tech)
            ]
            required_sheets = list(set(required_sheets))
            data = pd.read_excel(fn, sheet_name=required_sheets)
        else:
            data = pd.read_excel(fn, sheet_name=None)

        return {node: data}

    except Exception as e:
        raise Exception(
            f"Error reading PEMMDB data at {node} for climate year {cyear} and planning year {pyear}: {e}"
        )


def _process_thermal_capacities(
    node_tech_data: pd.DataFrame, node: str, cyear: int, pemmdb_tech: str
) -> pd.DataFrame:
    """
    Extract and clean thermal (conventionals & hydrogen) capacities.
    """
    # Extract capacities (p_nom)
    df = (
        node_tech_data.iloc[10:, [0, 1, 2]]
        .set_axis(["pemmdb_carrier", "pemmdb_type", "p_nom"], axis="columns")
        .dropna(how="all")
        .assign(
            pemmdb_carrier=lambda x: x.pemmdb_carrier.ffill(),
            pemmdb_type=lambda x: x.pemmdb_type.fillna(x.pemmdb_carrier),
            p_nom=lambda x: pd.to_numeric(x.p_nom, errors="coerce"),
            efficiency=1.0,  # capacities are in MWel and no efficiencies are given
            bus=node,
            country=node[:2],
            unit="MW",
        )
        .query("pemmdb_carrier == @pemmdb_tech")
        .reset_index(drop=True)
    )

    if df.empty:
        logger.info(
            f"No PEMMDB capacities match climate year {cyear} for '{pemmdb_tech}' at {node}."
        )
        return None

    return df


def _process_other_nonres_capacities(
    node_tech_data: pd.DataFrame, node: str, cyear: int, pemmdb_tech: str
) -> pd.DataFrame:
    """
    Extract and clean `Other Non-RES` capacities.
    """
    # Get raw data and filter for price band columns
    df = node_tech_data.iloc[6:16].reset_index(drop=True)

    df = (
        df.set_axis(df.iloc[0], axis="columns")
        .iloc[1:]
        .dropna(how="all")
        .filter(like="Price Band")
    )

    if df.empty:
        logger.info(
            f"No PEMMDB capacities available for '{pemmdb_tech}' and climate year {cyear} at node {node}."
        )
        return None

    column_names = [
        "p_nom",
        "units_count",
        "pemmdb_type",
        "purpose",
        "price",
        "efficiency",
        "co2_factor",
        "cyear_start",
        "cyear_end",
    ]

    # Extract data for given cyear
    df = (
        df.set_axis(column_names)
        .T.assign(
            pemmdb_carrier="Other Non-RES",
            bus=node,
            country=node[:2],
            unit="MW",
            pemmdb_type=lambda x: (
                x.pemmdb_type.str.split("/").str[1:].str.join("-")
                + "-"
                + x.purpose.astype(str)
                + "-"
                + x.price.astype(str)
                + "eur"
            ),
            cyear_start=lambda x: pd.to_numeric(x.cyear_start, errors="coerce"),
            cyear_end=lambda x: pd.to_numeric(x.cyear_end, errors="coerce"),
            p_nom=lambda x: pd.to_numeric(x.p_nom, errors="coerce"),
            units_count=lambda x: pd.to_numeric(x.units_count, errors="coerce"),
            price=lambda x: pd.to_numeric(x.price, errors="coerce"),
            efficiency=lambda x: pd.to_numeric(x.efficiency, errors="coerce"),
            co2_factor=lambda x: pd.to_numeric(x.co2_factor, errors="coerce"),
        )
        .query("cyear_start <= @cyear and cyear_end >= @cyear and p_nom > 0")
        .reset_index(drop=True)
    )

    if df.empty:
        logger.info(
            f"No PEMMDB capacity data matches climate year {cyear} for '{pemmdb_tech}' at {node}."
        )
        return None

    # Some datasets have duplicate pemmdb_tech price bands with same cyear, type, purpose and price but different capacities. Using first entry.
    if df.pemmdb_type.duplicated().any():
        logger.warning(
            f"Found duplicate '{pemmdb_tech}' price bands at {node} (cyear {cyear}) with same type, purpose, and price but different capacities. Using first entry."
        )

    df = df.groupby("pemmdb_type", as_index=False).first().reset_index(drop=True)

    return df


def _parse_index_parts(
    index: pd.Index, pemmdb_tech: str
) -> tuple[pd.Series, pd.Series]:
    """
    Parse index parts of PEMMDB renewable capacity sheets based on technology type.
    """
    if pemmdb_tech == "Hydro":
        pattern = r"^(.*?)\s*(?:\(.*?\))?\s*-\s*.*?capacity\s*\((.*?)\)"
    else:
        pattern = r"capacities\s+(.*?)\s+\((.*?)\)"
    extracted = index.str.extract(pattern).T.values

    return extracted[0], extracted[1]


def _process_res_capacities(
    node_tech_data: pd.DataFrame,
    node: str,
    cyear: int,
    pemmdb_tech: str,
) -> pd.DataFrame:
    """
    Extract and clean `RES` (Solar, Wind, Hydro) capacities.
    """
    # Clean data
    df = (
        node_tech_data.iloc[5:]
        .set_axis(["attributes", "p_nom"], axis="columns")
        .set_index("attributes")
        .rename_axis(None, axis=0)
        .dropna()
    )

    if df.empty:
        logger.info(
            f"No PEMMDB capacities match climate year {cyear} for '{pemmdb_tech}' at {node}."
        )
        return None

    # Infer type and unit from index
    types, units = _parse_index_parts(df.index, pemmdb_tech)

    # Extract data
    df = df.assign(
        bus=node,
        country=node[:2],
        pemmdb_carrier=pemmdb_tech,
        efficiency=1.0,  # no efficiencies given but capacity in MWel
        pemmdb_type=types,
        unit=units,
        element=lambda x: np.select(
            [
                x.index.str.contains("turbining", case=False, na=False),
                x.index.str.contains("pumping", case=False, na=False),
                x.index.str.contains("reservoir", case=False, na=False),
                x.index.str.contains("Storage capacities", case=False, na=False),
            ],
            ["Turbine", "Pump", "Reservoir", "Storage"],
            default="",
        ),
    )

    # Update type to include element information where needed
    df["pemmdb_type"] = np.where(
        df["element"] != "",
        df["pemmdb_type"] + " - " + df["element"],
        df["pemmdb_type"],
    )

    df = convert_units(df, value_col="p_nom").reset_index(drop=True)

    return df


def _process_other_res_capacities(
    node_tech_data: pd.DataFrame, node: str, cyear: int, pemmdb_tech: str
) -> pd.DataFrame:
    """
    Extract and clean `Other RES` capacities.
    """
    # Get raw data and extract capacity information
    df = (
        node_tech_data.iloc[[6], [1]]
        .set_axis(["p_nom"], axis=1)
        .reset_index(drop=True)
        .assign(
            bus=node,
            country=node[:2],
            pemmdb_carrier=pemmdb_tech,
            efficiency=1.0,  # no efficiencies given but capacity in MWel
            pemmdb_type="Small Biomass, Geothermal, Marine, Waste and Not Defined",
            unit="MW",
        )
        .dropna()
    )

    if df.empty:
        logger.info(
            f"No PEMMDB capacities available for '{pemmdb_tech}' and climate year {cyear} at node {node}."
        )
        return None

    return df


def _process_electrolyser_capacities(
    node_tech_data: pd.DataFrame, node: str, cyear: int, pemmdb_tech: str
) -> pd.DataFrame:
    """
    Extract and clean `Electrolyser` capacities.
    """
    # Extract data
    df = node_tech_data.iloc[7:].dropna(how="all", axis=0).dropna(how="all", axis=1)

    if df.empty:
        logger.info(
            f"No PEMMDB capacities available for '{pemmdb_tech}' and climate year {cyear} at node {node}."
        )
        return None

    column_names = [
        "attributes",
        "p_nom",
        "units_count",
        "efficiency",
        "h2_storage",
        "ramp_limit_up",
        "ramp_limit_down",
        "generation_reduction",
    ]

    df = (
        df.set_axis(column_names, axis=1)
        .set_index("attributes")
        .assign(
            pemmdb_carrier=pemmdb_tech,
            bus=node,
            country=node[:2],
            pemmdb_type="Onshore grid connected",
            unit="MW",
        )
        .reset_index(drop=True)
    )

    return df


def _process_battery_capacities(
    node_tech_data: pd.DataFrame, node: str, cyear: int, pemmdb_tech: str
) -> pd.DataFrame:
    """
    Extract and clean `Battery` capacities.
    """
    # Extract data
    df_raw = (
        node_tech_data.iloc[7:, 1:]
        .dropna(how="all", axis=0)
        .dropna(how="all", axis=1)
        .reset_index(drop=True)
    )

    if df_raw.empty:
        logger.info(
            f"No PEMMDB data available for '{pemmdb_tech}' and climate year {cyear} at node {node}."
        )
        return None

    column_names = [
        "p_nom_discharge",
        "p_nom_charge",
        "p_nom_store",
        "units_count",
        "efficiency",
        "ramp_limit_up",
        "ramp_limit_down",
    ]

    df_raw = df_raw.set_axis(column_names, axis=1)

    units = ["MW", "MW", "MWh"]
    types = ["Charge", "Discharge", "Store"]
    p_noms = pd.concat(
        [df_raw["p_nom_charge"], df_raw["p_nom_discharge"], df_raw["p_nom_store"]],
        axis=0,
    )

    df = pd.DataFrame(
        dict(
            p_nom=p_noms.values,
            efficiency=df_raw.efficiency[0],
            pemmdb_carrier=pemmdb_tech,
            bus=node,
            country=node[:2],
            pemmdb_type=types,
            unit=units,
        )
    )

    return df


def _process_dsr_capacities(
    node_tech_data: pd.DataFrame, node: str, cyear: int, pemmdb_tech: str
) -> pd.DataFrame:
    """
    Extract and clean `DSR` capacities.

    nb: Include only Market Response (user not willing to pay more than 'Activation price')
    nb: Activation price of -1 marks bands not to be modelled in MARKET models
    """
    # Get raw data
    df = node_tech_data.iloc[7:13, 2:]

    if df.empty:
        logger.info(
            f"No PEMMDB capacities available for '{pemmdb_tech}' and climate year {cyear} at node {node}."
        )
        return None

    column_names = [
        "p_nom",
        "units_count",
        "hours",
        "price",
        "cyear_start",
        "cyear_end",
    ]

    # Extract information and filter for given cyear
    df = (
        df.set_axis(column_names)
        .T.assign(
            pemmdb_carrier=pemmdb_tech,
            bus=node,
            country=node[:2],
            unit="MW",
            cyear_start=lambda x: pd.to_numeric(x.cyear_start, errors="coerce"),
            cyear_end=lambda x: pd.to_numeric(x.cyear_end, errors="coerce"),
            p_nom=lambda x: pd.to_numeric(x.p_nom, errors="coerce"),
            units_count=lambda x: pd.to_numeric(x.units_count, errors="coerce"),
            price=lambda x: pd.to_numeric(x.price, errors="coerce"),
            pemmdb_type=lambda x: x.hours.astype("str")
            + "h-"
            + x.price.astype("str")
            + "eur",
            efficiency=1.0,  # dummy value for efficiency
        )
        .query("cyear_start <= @cyear and cyear_end >= @cyear and p_nom > 0")
        .reset_index(drop=True)
    )

    if df.empty:
        logger.info(
            f"No PEMMDB capacity data matches climate year {cyear} for '{pemmdb_tech}' at {node}."
        )
        return None

    if df.pemmdb_type.duplicated().any():
        logger.warning(
            f"{node} has duplicated price bands for 'DSR' and cyear {cyear} with the same type (hours and price) but differing capacities. Keeping only first entry."
        )
        df = df.groupby("pemmdb_type", as_index=False).first()

    return df


def _process_thermal_profiles(
    node_tech_data: pd.DataFrame,
    node: str,
    pemmdb_tech: str,
    tyndp_scenario: str,
    pyear_i: int,
    sns: pd.DatetimeIndex,
    index_year: pd.DatetimeIndex,
) -> pd.DataFrame:
    """
    Extract and clean thermal (conventionals & hydrogen) profiles.
    """
    # Extract data
    must_runs = (
        node_tech_data.iloc[10:, [0, 1, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]]
        .set_axis(
            [
                "pemmdb_carrier",
                "pemmdb_type",
                "unit",
                "Jan",
                "Feb",
                "Mar",
                "Apr",
                "May",
                "Jun",
                "Jul",
                "Aug",
                "Sep",
                "Oct",
                "Nov",
                "Dec",
            ],
            axis="columns",
        )
        .dropna(how="all")
    )

    # Deal with Nuclear and Light Oil entries as they do not have a PEMMDB type
    nuclear_lightoil_i = must_runs.query(
        "pemmdb_carrier in ['Nuclear', 'Light oil']"
    ).index
    must_runs.loc[nuclear_lightoil_i, "pemmdb_type"] = must_runs.loc[
        nuclear_lightoil_i, "pemmdb_carrier"
    ]
    must_runs = (
        must_runs.assign(
            pemmdb_carrier=lambda df: df.pemmdb_carrier.ffill(),
            pemmdb_type=lambda df: df.pemmdb_type.ffill(),
        )
        .query("unit == '% of installed capacity' and pemmdb_carrier == @pemmdb_tech")
        .drop(columns=["unit"])
        .set_index(["pemmdb_type"])
        .drop(columns=["pemmdb_carrier"])
        .div(1e2)  # percentage conversion
    )

    # Transform to monthly lookup by inverting df
    must_runs = must_runs.T

    # Map to hourly Datetime index
    df = pd.DataFrame({"month": index_year.month_name().str[:3]}, index=index_year)
    df = df.join(must_runs, on="month")
    df.drop(columns="month", inplace=True)

    # Flatten and turn into xarray dataset
    profiles = (
        df.loc[sns]
        .melt(var_name="pemmdb_type", value_name="p_min_pu", ignore_index=False)
        .rename_axis("time", axis="index")
        .assign(
            pemmdb_carrier=pemmdb_tech, bus=node, p_max_pu=1.0
        )  # also set p_max_pu with default value of 1.0
        .set_index(["bus", "pemmdb_carrier", "pemmdb_type"], append=True)
        .sort_index()  # sort index for more efficient indexing
    )

    # Remove must-runs for DE and GA scenarios after 2030 (per TYNDP 2024 Methodology, p.37)
    if tyndp_scenario != "NT" and pyear_i > 2030:
        profiles.loc[:, "p_min_pu"] = 0.0

    return profiles


def _process_other_res_profiles(
    node_tech_data: pd.DataFrame,
    node: str,
    pemmdb_tech: str,
    pyear: int,
    sns: pd.DatetimeIndex,
    index_year: pd.DatetimeIndex,
) -> pd.DataFrame:
    """
    Extract and clean `Other RES` profiles.
    """
    # Extract data
    df = node_tech_data.iloc[6:, :3]

    capacity = np.float64(df.iloc[0, 1])

    if np.isnan(capacity):
        logger.warning(
            f"No 'Other RES' capacity found for {node} in {pyear}, hence no must-run profile available."
        )
        return None

    profiles = (
        df.iloc[3:, 2]
        .to_frame()
        .set_axis(["p_min"], axis="columns")
        .assign(
            p_min_pu=lambda df: pd.to_numeric(df.p_min / capacity)
            if capacity > 0
            else 0.0,
            p_max_pu=1.0,  # also set p_max_pu with default value of 1.0
            time=index_year,
            bus=node,
            pemmdb_carrier=pemmdb_tech,
            pemmdb_type="Small Biomass, Geothermal, Marine, Waste and Not Defined",
        )
        .loc[lambda x: x["time"].isin(sns)]
        .set_index(["time", "bus", "pemmdb_carrier", "pemmdb_type"])[
            ["p_min_pu", "p_max_pu"]
        ]
    )

    return profiles


def _process_other_nonres_profiles(
    node_tech_data: pd.DataFrame,
    node: str,
    pemmdb_tech: str,
    cyear: int,
    sns: pd.DatetimeIndex,
    index_year: pd.DatetimeIndex,
) -> pd.DataFrame:
    """
    Extract and clean `Other Non-RES` profiles.
    """
    # Extract data
    df = node_tech_data.iloc[7:, 1:]
    df = (
        df.set_index(df.columns[0])
        .rename(
            index={
                "Installed capacity \n(MW)": "p_nom",
                "PEMMDB type(s)": "pemmdb_type",
                "Purpose": "purpose",
                "Avg. Market Offer Price (€/MWh)": "price",
                "Avg. Efficiency Ratio": "efficiency",
                "Start climate year": "cyear_start",
                "End climate year": "cyear_end",
            }
        )
        .rename_axis(None, axis=0)
    )

    # Create mask to filter for given climate year
    cyear_start = pd.to_numeric(df.loc["cyear_start", :], errors="coerce")
    cyear_end = pd.to_numeric(df.loc["cyear_end", :], errors="coerce")
    cap = pd.to_numeric(df.loc["p_nom", :], errors="coerce")
    mask = (cyear_start <= cyear) & (cyear <= cyear_end) & (cap > 0)

    if not mask.any():
        logger.warning(
            f"No PEMMDB profiles available for '{pemmdb_tech}' and climate year {cyear} at node {node}."
        )
        return None

    # Filter for climate year
    df = df.loc[:, mask]

    # Extract capacity and plant type
    type = (
        df.loc["pemmdb_type", :].str.split("/").str[1:].str.join("-")
        + "-"
        + df.loc["purpose", :].astype(str)
        + "-"
        + df.loc["price", :].astype(str)
        + "eur"
    )

    df_long = (
        df.iloc[10:]
        .reset_index(drop=True)
        .div(cap.loc[mask], axis=1)
        .set_axis(type, axis="columns")
        .assign(
            time=index_year,
            bus=node,
            pemmdb_carrier=pemmdb_tech,
        )
        .loc[lambda x: x["time"].isin(sns)]
        .set_index(["time", "bus", "pemmdb_carrier"])
        .melt(var_name="pemmdb_type", value_name="p_max_pu", ignore_index=False)
        .set_index("pemmdb_type", append=True)
        .assign(p_min_pu=0.0)  # also set p_min_pu with default value of 0.0
    )

    if df_long.index.duplicated().any():
        logger.warning(
            f"{node} has duplicated price bands for '{pemmdb_tech}' and cyear {cyear} with the same type (type, purpose, price) but differing capacities. Keeping only first entry."
        )

    # Some datasets have duplicate price bands with same cyear, hours and price but different capacities. We keep the first entry only
    profiles = df_long.groupby(df_long.index.names).first()

    return profiles


def _process_dsr_profiles(
    node_tech_data: pd.DataFrame,
    node: str,
    pemmdb_tech: str,
    cyear: int,
    sns: pd.DatetimeIndex,
    index_year: pd.DatetimeIndex,
) -> pd.DataFrame:
    """
    Extract and clean `DSR` profiles.
    """
    # Extract data
    df = node_tech_data.iloc[7:, 1:]
    df = df.set_index(df.columns[0]).rename(
        index={
            "Capacity": "p_nom",
            "Units": "units_count",
            "Hours": "hours",
            "Price": "price",
            "Climate year start": "cyear_start",
            "Climate year end": "cyear_end",
        }
    )

    # Create mask to filter for given climate year and for capacity > 0
    cyear_start = pd.to_numeric(df.loc["cyear_start", :], errors="coerce")
    cyear_end = pd.to_numeric(df.loc["cyear_end", :], errors="coerce")
    cap = pd.to_numeric(df.loc["p_nom", :], errors="coerce")
    mask = (cyear_start <= cyear) & (cyear <= cyear_end) & (cap > 0)

    if not mask.any():
        logger.warning(
            f"No PEMMDB data available for '{pemmdb_tech}' and climate year {cyear} at node {node}."
        )
        return None

    # Filter for climate year
    df = df.loc[:, mask]

    # Extract price band type information
    type = (
        df.loc["hours", :].astype("str")
        + "h-"
        + df.loc["price", :].astype("str")
        + "eur"
    )

    df_long = (
        df.iloc[7:]
        .reset_index(drop=True)
        .div(cap.loc[mask], axis=1)
        .set_axis(type, axis="columns")
        .assign(
            time=index_year,
            bus=node,
            pemmdb_carrier=pemmdb_tech,
        )
        .loc[lambda x: x["time"].isin(sns)]
        .set_index(["time", "bus", "pemmdb_carrier"])
        .melt(var_name="pemmdb_type", value_name="p_max_pu", ignore_index=False)
        .set_index("pemmdb_type", append=True)
        .assign(p_min_pu=0.0)  # also set p_min_pu with default value of 0.0
    )

    if df_long.index.duplicated().any():
        logger.warning(
            f"{node} has duplicated price bands for 'DSR' and cyear {cyear} with the same type (hours and price) but differing capacities. Keeping only first entry."
        )

    # Some datasets have duplicate pemmdb_tech price bands with same cyear, type, purpose and price but different capacities. Using first entry.
    profiles = df_long.groupby(df_long.index.names).first()

    return profiles


def process_pemmdb_capacities(
    node_tech_data: pd.DataFrame,
    node: str,
    pemmdb_tech: str,
    cyear: int,
    pyear: int,
    carrier_mapping_fn: str,
) -> pd.DataFrame:
    """
    Read and clean capacities from PEMMDB for a given technology, planning and climate year.

    Parameters
    ----------
    node_tech_data : pd.DataFrame
        Dataframe containing the PEMMDB capacity data sheet for the given node and technology.
    node : str
        Node name to read data for.
    pemmdb_tech : str
        PEMMDB technology to read data for.
    cyear : int
        Climate year to read data for.
    pyear : int
        Planning year used for data retrieval (fallback year if pyear_i not available).
    carrier_mapping_fn : str
        Path to file with mapping from external carriers to available tyndp_carrier names.

    Returns
    -------
    pd.DataFrame
        Dataframe containing NT capacities (p_nom) in long format for the given PEMMDB technology and node.
    """
    try:
        # Conventionals & Hydrogen
        if pemmdb_tech in CONVENTIONALS or pemmdb_tech == "Hydrogen":
            capacities = _process_thermal_capacities(
                node_tech_data, node, cyear, pemmdb_tech
            )

        # Other Non-RES
        elif pemmdb_tech == "Other Non-RES":
            capacities = _process_other_nonres_capacities(
                node_tech_data, node, cyear, pemmdb_tech
            )

        # Renewables (Solar, Wind, Hydro)
        elif pemmdb_tech in RENEWABLES:
            capacities = _process_res_capacities(
                node_tech_data, node, cyear, pemmdb_tech
            )

        # Other RES
        elif pemmdb_tech == "Other RES":
            capacities = _process_other_res_capacities(
                node_tech_data, node, cyear, pemmdb_tech
            )

        # DSR
        elif pemmdb_tech == "DSR":
            capacities = _process_dsr_capacities(
                node_tech_data, node, cyear, pemmdb_tech
            )

        # Battery
        elif pemmdb_tech == "Battery":
            capacities = _process_battery_capacities(
                node_tech_data, node, cyear, pemmdb_tech
            )

        # Electrolyser
        elif pemmdb_tech == "Electrolyser":
            capacities = _process_electrolyser_capacities(
                node_tech_data, node, cyear, pemmdb_tech
            )

        else:
            return None

        if capacities is None:
            return None

        # Separate energy and power capacities and select needed columns
        capacities = (
            capacities.assign(
                e_nom=lambda x: np.where(x.unit.str.contains("h"), x.p_nom, 0.0),
                p_nom=lambda x: np.where(x.unit.str.contains("h"), 0.0, x.p_nom),
            )[
                [
                    "pemmdb_carrier",
                    "bus",
                    "pemmdb_type",
                    "p_nom",
                    "e_nom",
                    "efficiency",
                    "country",
                    "unit",
                ]
            ]  # select and order relevant columns
        ).reset_index(drop=True)

        # Map pemmdb_carrier and pemmdb_type to TYNDP technology names
        capacities = map_tyndp_carrier_names(
            capacities,
            carrier_mapping_fn,
            ["pemmdb_carrier", "pemmdb_type"],
            drop_merge_columns=True,
        )

        return capacities

    except Exception as e:
        raise Exception(
            f"Error while processing capacities for {pemmdb_tech} at {node} for climate year {cyear} and planning year {pyear}: {e}"
        )


def process_pemmdb_profiles(
    node_tech_data: pd.DataFrame,
    node: str,
    pemmdb_tech: str,
    tyndp_scenario: str,
    cyear: int,
    pyear: int,
    pyear_i: int,
    sns: pd.DatetimeIndex,
    index_year: pd.DatetimeIndex,
    carrier_mapping_fn: str,
) -> pd.DataFrame:
    """
    Read and clean must run obligations (p_min_pu) and availability (p_max_pu) profiles
    from PEMMDB for a given technology, planning and climate year.

    Parameters
    ----------
    node_tech_data : pd.DataFrame
        Dataframe containing the PEMMDB profile data sheet for the given node and technology.
    node : str
        Node name to read data for.
    pemmdb_tech : str
        PEMMDB technology to read data for.
    tyndp_scenario : str
        TYNDP scenario to read data for.
    cyear : int
        Climate year to read data for.
    pyear : int
        Planning year used for data retrieval (fallback year if pyear_i not available).
    pyear_i : int
        Original planning year.
    sns : pd.DatetimeIndex
        Modelled snapshots.
    index_year : pd.DatetimeIndex
        Hourly Datetime index for a full given cyear.
    carrier_mapping_fn : str
        Path to file with mapping from external carriers to available tyndp_carrier names

    Returns
    -------
    profiles : pd.DataFrame
        Dataframe containing must run obligations (p_min_pu) and availability (p_max_pu) profiles.
    """
    try:
        # Conventionals & Hydrogen
        if pemmdb_tech in CONVENTIONALS or pemmdb_tech == "Hydrogen":
            profiles = _process_thermal_profiles(
                node_tech_data,
                node,
                pemmdb_tech,
                tyndp_scenario,
                pyear_i,
                sns,
                index_year,
            )

        # Other RES
        elif pemmdb_tech == "Other RES":
            profiles = _process_other_res_profiles(
                node_tech_data, node, pemmdb_tech, pyear, sns, index_year
            )

        # Other Non-RES
        elif pemmdb_tech == "Other Non-RES":
            profiles = _process_other_nonres_profiles(
                node_tech_data, node, pemmdb_tech, cyear, sns, index_year
            )

        # DSR
        elif pemmdb_tech == "DSR":
            profiles = _process_dsr_profiles(
                node_tech_data, node, pemmdb_tech, cyear, sns, index_year
            )

        else:
            return None

        if profiles is None:
            return None

        # Map PEMMDB carrier names to TYNDP technologies
        profiles = map_tyndp_carrier_names(
            profiles.reset_index(),
            carrier_mapping_fn,
            ["pemmdb_carrier", "pemmdb_type"],
            drop_merge_columns=True,
        ).set_index(["time", "bus", "carrier", "index_carrier"])

        return profiles

    except Exception as e:
        raise Exception(
            f"Error reading PEMMDB profiles for '{pemmdb_tech}' at {node} for climate year {cyear} and planning year {pyear}: {e}"
        )


def process_pemmdb_data(
    category: str,
    node_tech: tuple[str, str],
    pemmdb_data: dict[str, dict[str, pd.DataFrame]],
    cyear: int,
    pyear: int,
    pyear_i: int,
    tyndp_scenario: str,
    sns: pd.DatetimeIndex,
    index_year: pd.DatetimeIndex,
    carrier_mapping_fn: str,
) -> pd.DataFrame:
    """
    Read and clean either capacities or must run obligations (p_min_pu) and availability (p_max_pu) profiles
    from PEMMDB for a given technology, planning and climate year.

    Parameters
    ----------
    category : str
        Category to read data for. Can be either 'capacities' or 'profiles'.
    node_tech : tuple[str, str]
        Tuple with node and technology to process for.
    pemmdb_data : dict[str, dict[str, pd.DataFrame]]
        Dictionary containing all PEMMDB data for all nodes and technologies.
    cyear : int
        Climate year to read data for.
    pyear : int
        Planning year used for data retrieval (fallback year if pyear_i not available).
    pyear_i : int
        Original planning year.
    tyndp_scenario : str
        TYNDP scenario to read data for.
    sns : pd.DatetimeIndex
        Modelled snapshots.
    index_year : pd.DatetimeIndex
        Hourly Datetime index for a full given cyear.
    carrier_mapping_fn : str
        Path to file with mapping from external carriers to available tyndp_carrier names.

    Returns
    -------
    pd.DataFrame
        Pandas dataframe containing capacities or must run obligations (p_min_pu) and availability (p_max_pu) profiles respectively.
    """
    # Extract node and tech information
    node, pemmdb_tech = node_tech

    # Extract PEMMDB data for corresponding node and tech
    node_tech_data = pemmdb_data.get(node, {}).get(
        PEMMDB_SHEET_MAPPING.get(pemmdb_tech, pemmdb_tech), None
    )
    if node_tech_data is None:
        return None

    if category == "capacities":
        data = process_pemmdb_capacities(
            node_tech_data,
            node,
            pemmdb_tech,
            cyear,
            pyear,
            carrier_mapping_fn,
        )
    elif category == "profiles":
        data = process_pemmdb_profiles(
            node_tech_data,
            node,
            pemmdb_tech,
            tyndp_scenario,
            cyear,
            pyear,
            pyear_i,
            sns,
            index_year,
            carrier_mapping_fn,
        )
    else:
        raise Exception(
            f"Unknown element for PEMMDB data: {category}. Please choose between 'capacities' and 'profiles'."
        )

    return data


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_pemmdb_data",
            clusters="all",
            planning_horizons=2030,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameter
    pemmdb_techs = [tech.replace("_", " ") for tech in snakemake.params.pemmdb_techs]
    nodes = pd.read_csv(snakemake.input.busmap, index_col=0).index
    pemmdb_dir = snakemake.input.pemmdb_dir
    tyndp_scenario = snakemake.params.tyndp_scenario
    carrier_mapping_fn = snakemake.input.carrier_mapping

    # Climate year from snapshots
    sns = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)
    cyear = sns[0].year
    index_year = pd.date_range(
        start=f"{cyear}-01-01",
        periods=8760,
        freq="h",
    )

    # Only climate years 1995, 2008 and 2009 are available for all technologies and countries
    if cyear not in [1995, 2008, 2009]:
        logger.warning(
            f"Snapshot year {cyear} doesn't match available TYNDP data. Falling back to 2009."
        )
        cyear = 2009

    # Planning year
    pyear_i = int(snakemake.wildcards.planning_horizons)
    pyear = safe_pyear(
        pyear_i,
        available_years=snakemake.params.available_years,
        source="PEMMDB",
    )

    # Load all PEMMDB data
    tqdm_kwargs_read = {
        "ascii": False,
        "unit": " nodes",
        "total": len(nodes),
        "desc": "Loading PEMMDB data...",
    }

    func_read = partial(
        read_pemmdb_data,
        pemmdb_dir=pemmdb_dir,
        cyear=cyear,
        pyear=pyear,
        required_techs=pemmdb_techs,
    )

    with mp.Pool(processes=snakemake.threads) as pool1:
        pemmdb_data_list = [
            data
            for data in tqdm(pool1.imap(func_read, nodes), **tqdm_kwargs_read)
            if data is not None
        ]

    pemmdb_data = {node: data for d in pemmdb_data_list for node, data in d.items()}

    ####################
    # Process capacities
    ####################

    node_techs = list(product(nodes, pemmdb_techs))

    tqdm_kwargs_caps = {
        "ascii": False,
        "unit": " node-tech combinations",
        "total": len(node_techs),
        "desc": "Processing PEMMDB capacities...",
    }

    with logging_redirect_tqdm():
        pemmdb_capacities = [
            caps
            for caps in (
                process_pemmdb_data(
                    "capacities",
                    node_tech=node_tech,
                    pemmdb_data=pemmdb_data,
                    cyear=cyear,
                    pyear=pyear,
                    pyear_i=pyear_i,
                    tyndp_scenario=tyndp_scenario,
                    sns=sns,
                    index_year=index_year,
                    carrier_mapping_fn=carrier_mapping_fn,
                )
                for node_tech in tqdm(node_techs, **tqdm_kwargs_caps)
            )
            if caps is not None
        ]

    if not pemmdb_capacities:
        logger.warning(
            f"No PEMMDB capacities available for climate year {cyear} and planning year {pyear}. "
            f"Please specify different technologies, climate year or planning year."
        )
        # Save empty file
        pd.DataFrame().to_csv(snakemake.output.pemmdb_capacities)
        sys.exit(0)

    # Otherwise concat capacities into one dataframe and save to csv
    pemmdb_capacities_df = pd.concat(pemmdb_capacities, axis=0)
    pemmdb_capacities_df.to_csv(snakemake.output.pemmdb_capacities, index=False)

    ##################
    # Process profiles
    ##################
    #
    # # Load and prep data
    tqdm_kwargs_profiles = {
        "ascii": False,
        "unit": " node-tech combinations",
        "total": len(node_techs),
        "desc": "Processing PEMMDB profiles...",
    }
    with logging_redirect_tqdm():
        pemmdb_profiles = [
            profiles
            for profiles in (
                process_pemmdb_data(
                    "profiles",
                    node_tech=node_tech,
                    pemmdb_data=pemmdb_data,
                    cyear=cyear,
                    pyear=pyear,
                    pyear_i=pyear_i,
                    tyndp_scenario=tyndp_scenario,
                    sns=sns,
                    index_year=index_year,
                    carrier_mapping_fn=carrier_mapping_fn,
                )
                for node_tech in tqdm(node_techs, **tqdm_kwargs_profiles)
            )
            if profiles is not None
        ]

    if not pemmdb_profiles:
        logger.warning(
            f"No PEMMDB profiles available for climate year {cyear} and planning year {pyear}. "
            f"Please specify different technologies, climate year or planning year."
        )
        # Save empty dataset
        xr.Dataset().to_netcdf(snakemake.output.pemmdb_profiles)
        sys.exit(0)

    # Otherwise merge PEMMDB profiles into one pd.DataFrame, convert to xarray dataset and save
    ds = pd.concat(pemmdb_profiles, axis=0).to_xarray()
    ds.to_netcdf(snakemake.output.pemmdb_profiles)
