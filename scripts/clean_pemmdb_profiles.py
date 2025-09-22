# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT
"""
Loads and cleans the available PEMMDB v2.4 must-run and availability profiles from TYNDP data bundle for a given

* climate year,
* planning horizon,
* technology.

Available technologies are:

* Conventionals
    * Nuclear
    * Hard coal
    * Lignite
    * Gas
    * Light oil
    * Heavy oil
    * Oil shale
    * Other Non-RES
* Renewables
    * Wind
    * Solar
    * Hydro
    * Other RES
* Hydrogen
    * Hydrogen fuel cell
    * Hydrogen CCGT
* Reserves
* DSR
* Battery
* Electrolyser


Outputs
-------
Cleaned netcdf file with must run obligations (p_min_pu) and availability (p_max_pu) profiles.
"""

import logging
import multiprocessing as mp
import sys
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm

from scripts._helpers import (
    configure_logging,
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


def _read_thermal_profiles(
    fn: Path,
    pemmdb_tech: str,
    tyndp_scenario: str,
    pyear_i: int,
    node: str,
    sns: pd.DatetimeIndex,
    index_year: pd.DatetimeIndex,
) -> xr.Dataset:
    """Read and clean thermal (conventionals & hydrogen) profiles."""

    # read data
    must_runs = (
        pd.read_excel(
            fn,
            sheet_name="Thermal",
            skiprows=7,
            usecols=[0, 1, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18],
            names=[
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
        )
        .drop([0, 1, 2])
        .dropna(how="all")
    )
    # deal with Nuclear and Light Oil entries as they do not have a pemmdb type
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

    # map to hourly Datetime index
    df = pd.DataFrame({"month": index_year.month_name().str[:3]}, index=index_year)
    for type in must_runs.columns:
        df[type] = df["month"].map(must_runs[type])
    df.drop(columns="month", inplace=True)

    # flatten and turn into xarray dataset
    profiles = (
        df.reset_index(names="time")
        .melt(id_vars="time", var_name="pemmdb_type", value_name="p_min_pu")
        .assign(
            pemmdb_carrier=pemmdb_tech, bus=node, p_max_pu=1.0
        )  # also set p_max_pu with default value of 1.0
        .set_index(["time", "bus", "pemmdb_carrier", "pemmdb_type"])
        .to_xarray()
        .sel(time=sns)  # filter for snapshots only
    )

    # Remove must-runs for DE and GA after 2030 as to 2024 TYNDP Methodology report, p.37
    if tyndp_scenario != "NT" and pyear_i > 2030:
        profiles["p_min_pu"] = xr.full_like(profiles["p_min_pu"], 0.0)

    return profiles


def _read_other_res_profiles(
    fn: Path,
    pemmdb_tech: str,
    node: str,
    pyear: int,
    sns: pd.DatetimeIndex,
    index_year: pd.DatetimeIndex,
) -> xr.Dataset:
    """Read and clean `Other RES` profiles."""

    # read data
    df = pd.read_excel(fn, sheet_name="Other RES", skiprows=6, usecols=[0, 1, 2])

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
        .set_index(["time", "bus", "pemmdb_carrier", "pemmdb_type"])[
            ["p_min_pu", "p_max_pu"]
        ]
        .to_xarray()
        .sel(time=sns)  # filter for snapshots only
    )

    return profiles


def _read_other_nonres_profiles(
    fn: Path,
    pemmdb_tech: str,
    cyear: int,
    node: str,
    sns: pd.DatetimeIndex,
    index_year: pd.DatetimeIndex,
) -> xr.Dataset:
    """
    Read and clean `Other Non-RES` profiles.
    """

    # read data
    df = pd.read_excel(fn, sheet_name="Other Non-RES", skiprows=7, index_col=1).rename(
        index={
            "Installed capacity \n(MW)": "p_nom",
            "PEMMDB type(s)": "pemmdb_type",
            "Start climate year": "cyear_start",
            "End climate year": "cyear_end",
        }
    )

    # Create mask to filter for given climate year
    cyear_start = pd.to_numeric(df.loc["cyear_start", :], errors="coerce")
    cyear_end = pd.to_numeric(df.loc["cyear_end", :], errors="coerce")
    mask = (cyear_start <= cyear) & (cyear <= cyear_end)

    if not mask.any():
        logger.warning(
            f"No PEMMDB data available for '{pemmdb_tech}' and climate year {cyear} at node {node}."
        )
        return None

    # Filter for climate year and rename single column
    df = df.loc[:, mask]
    if df.shape[1] != 1:
        logger.warning(
            f"Expected max 1 column to match climate year {cyear} for '{pemmdb_tech}' at {node}, got {df.shape[1]}. Using first entry only."
        )
        df = df.iloc[:, :1]  # Take only first column

    df = df.set_axis(["p_max"], axis="columns")

    # Extract capacity and plant type
    cap = pd.to_numeric(df.loc["p_nom", "p_max"], errors="coerce")
    type = " - ".join(df.loc["pemmdb_type", "p_max"].split("/")[1:])

    if pd.isna(cap) or cap <= 0:
        logger.warning(f"Invalid capacity data for '{pemmdb_tech}' at {node}")
        return None

    profiles = (
        df.iloc[10:]
        .reset_index(drop=True)
        .assign(
            p_nom=cap,
            p_max_pu=lambda x: x.p_max / cap,  # calculate per unit availability
            p_min_pu=0.0,  # also set p_min_pu with default value of 0.0
            time=index_year,
            bus=node,
            pemmdb_carrier=pemmdb_tech,
            pemmdb_type=type,
        )
        .set_index(["time", "bus", "pemmdb_carrier", "pemmdb_type"])[
            ["p_min_pu", "p_max_pu"]
        ]
        .to_xarray()
        .sel(time=sns)
    )

    return profiles


def _read_dsr_profiles(
    fn: Path,
    pemmdb_tech: str,
    cyear: int,
    node: str,
    sns: pd.DatetimeIndex,
    index_year: pd.DatetimeIndex,
) -> xr.Dataset:
    """
    Read and clean `DSR` profiles.
    """

    # read data
    df = pd.read_excel(fn, sheet_name="DSR", skiprows=7, index_col=1).rename(
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

    # Filter for climate year and rename single column
    df = df.loc[:, mask]

    # extract price band type information
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
        .set_index(["time", "bus", "pemmdb_carrier"])
        .melt(var_name="pemmdb_type", value_name="p_max_pu", ignore_index=False)
        .set_index("pemmdb_type", append=True)
        .assign(p_min_pu=0.0)  # also set p_min_pu with default value of 0.0
    )

    if df_long.index.duplicated().any():
        logger.warning(
            f"{node} has duplicated price bands for 'DSR' and cyear {cyear} with the same type (hours and price) but differing capacities. Keeping only first entry."
        )

    # some datasets have duplicate price bands with same cyear, hours and price but different capacities. We keep the first entry only
    profiles = df_long.groupby(level=[0, 1, 2, 3]).first().to_xarray().sel(time=sns)

    return profiles


def read_pemmdb_profiles(
    node: str,
    pemmdb_dir: str,
    tyndp_scenario: str,
    cyear: int,
    pyear: int,
    pyear_i: int,
    pemmdb_tech: str,
    sns: pd.DatetimeIndex,
    index_year: pd.DatetimeIndex,
) -> xr.Dataset:
    """
    Reads and cleans must run obligations (p_min_pu) and availability (p_max_pu) profiles
    from PEMMDB for a given technology, planning and climate year.

    Parameters
    ----------
    node : str
        Node name to read data for.
    pemmdb_dir : str
        Path to directory containing PEMMDB data.
    tyndp_scenario : str
        TYNDP scenario to read data for.
    cyear : int
        Climate year to read data for.
    pyear : int
        Planning year to read data for. Can be fallback year to available data.
    pyear_i : int
        Original planning year.
    pemmdb_tech : str
        PEMMDB technology to read data for.
    sns : pd.DatetimeIndex
        Modelled snapshots.
    index_year : pd.DatetimeIndex
        Hourly Datetime index for a full given cyear.

    Returns
    -------
    profiles : xr.Dataset
        xarray dataset containing must run obligations (p_min_pu) and availability (p_max_pu) profiles.
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
        # Conventionals & Hydrogen
        if pemmdb_tech in CONVENTIONALS or pemmdb_tech == "Hydrogen":
            return _read_thermal_profiles(
                fn, pemmdb_tech, tyndp_scenario, pyear_i, node, sns, index_year
            )

        # Other RES
        elif pemmdb_tech == "Other RES":
            return _read_other_res_profiles(
                fn, pemmdb_tech, node, pyear, sns, index_year
            )

        # Other Non-RES
        elif pemmdb_tech == "Other Non-RES":
            return _read_other_nonres_profiles(
                fn, pemmdb_tech, cyear, node, sns, index_year
            )

        # DSR
        elif pemmdb_tech == "DSR":
            return _read_dsr_profiles(fn, pemmdb_tech, cyear, node, sns, index_year)

        else:
            return None

    except Exception as e:
        raise Exception(
            f"Error reading PEMMDB profiles for '{pemmdb_tech}' at {node} for climate year {cyear} and planning year {pyear}: {e}"
        )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_pemmdb_profiles",
            clusters="all",
            planning_horizons=2030,
            tech="Gas",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Climate year from snapshots
    sns = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)
    cyear = sns[0].year
    index_year = pd.date_range(
        start=f"{cyear}-01-01",
        periods=8760,  # 53 weeks
        freq="h",
    )
    # Only climate years 1995, 2008 and 2009 are available for all technologies and countries
    if cyear not in [1995, 2008, 2009]:
        logger.warning(
            "Snapshot year doesn't match available TYNDP data. Falling back to 2009."
        )
        cyear = 2009

    # Planning year
    pyear_i = int(snakemake.wildcards.planning_horizons)
    pyear = safe_pyear(
        pyear_i,
        available_years=snakemake.params.available_years,
        source="PEMMDB",
    )

    # Parameters
    onshore_buses = pd.read_csv(snakemake.input.busmap, index_col=0)
    nodes = onshore_buses.index
    pemmdb_dir = snakemake.input.pemmdb_dir
    tech = str(snakemake.wildcards.tech)
    tyndp_scenario = snakemake.params.tyndp_scenario
    carrier_mapping_df = (
        pd.read_csv(snakemake.input.carrier_mapping)[
            ["pemmdb_carrier", "pemmdb_type", "open_tyndp_carrier", "open_tyndp_index"]
        ]
    ).dropna()

    # Load and prep data
    tqdm_kwargs = {
        "ascii": False,
        "unit": " nodes",
        "total": len(nodes),
        "desc": f"Loading PEMMDB profiles for {tech}",
    }

    func = partial(
        read_pemmdb_profiles,
        pemmdb_dir=pemmdb_dir,
        tyndp_scenario=tyndp_scenario,
        cyear=cyear,
        pyear=pyear,
        pyear_i=pyear_i,
        pemmdb_tech=tech,
        sns=sns,
        index_year=index_year,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        profiles = [
            profile
            for profile in tqdm(pool.imap(func, nodes), **tqdm_kwargs)
            if profile is not None
        ]

    if not profiles:
        logger.warning(
            f"No PEMMDB profiles available for '{tech}' with climate year {cyear} and planning year {pyear}. "
            f"Please specify different technology, climate year or planning year."
        )
        # save empty dataset
        ds = xr.Dataset()
        ds.to_netcdf(snakemake.output.pemmdb_profiles)
        sys.exit(0)

    # otherwise merge pemmdb profiles into one xarray dataset and save
    ds = xr.merge(profiles)
    # map pemmdb carrier names to tyndp technologies
    ds = (
        map_tyndp_carrier_names(
            ds.to_dataframe().reset_index(),
            carrier_mapping_df,
            ["pemmdb_carrier", "pemmdb_type"],
        )
        .set_index(["time", "bus", "carrier", "index_carrier"])
        .to_xarray()
    )
    ds.to_netcdf(snakemake.output.pemmdb_profiles)
