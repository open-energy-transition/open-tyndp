# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT
"""
Loads and cleans the available PEMMDB v2.4 data from TYNDP data bundle for a given

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
    * Other non-RES
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
Cleaned csv file with NT capacities (p_nom) in long format and netcdf file with must run obligations (p_min_pu) profiles.
"""

import logging
import multiprocessing as mp
import os
from functools import partial
from pathlib import Path

import pandas as pd
import xarray as xr
from tqdm import tqdm

from scripts._helpers import (
    configure_logging,
    get_snapshots,
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
    "Other non-RES",
]


def read_pemmdb_capacities(
    node: str,
    pemmdb_dir: str,
    cyear: str,
    pyear: str,
    pemmdb_tech: str,
    sns: pd.DatetimeIndex,
) -> pd.Series:
    fn = Path(
        pemmdb_dir,
        str(pyear),
        f"PEMMDB_{node.replace('GB', 'UK')}_NationalTrends_{pyear}.xlsx",
    )

    if not os.path.isfile(fn):
        return None

    # Conventionals & Hydrogen
    if pemmdb_tech in CONVENTIONALS or pemmdb_tech == "Hydrogen":
        # read capacities (p_nom)
        df = (
            pd.read_excel(
                fn,
                sheet_name="Thermal",
                skiprows=7,  # first seven rows are metadata
                usecols=[0, 1, 2],
                names=["carrier", "type", "p_nom"],
            )
            .drop([0, 1, 2])
            .dropna(how="all")
            .assign(
                **{
                    "carrier": lambda df: df["carrier"].ffill(),
                    "type": lambda df: df["type"].fillna(df.carrier),
                    "bus": node,
                    "country": node[:2],
                }
            )
            .query("carrier == @pemmdb_tech")
        )

    # Wind
    elif pemmdb_tech == "Wind":
        df = pd.read_excel(fn, sheet_name="Wind")

    # Solar
    elif pemmdb_tech == "Solar":
        df = pd.read_excel(
            fn,
            sheet_name="Solar",
            # skiprows=1,
            # usecols=lambda name: name == "Day"
            # or name == "Week"
            # or name == "ShortName"
            # or name == "Variable"
            # or name == int(cyear),
            # sheet_name=f"{hydro_tech} - Year Dependent",
        )

    # Hydro
    elif pemmdb_tech == "hydro":
        df = pd.read_excel(
            fn,
            sheet_name="Hydro",
            # skiprows=1,
            # usecols=lambda name: name == "Day"
            # or name == "Week"
            # or name == "ShortName"
            # or name == "Variable"
            # or name == int(cyear),
            # sheet_name=f"{hydro_tech} - Year Dependent",
        )

    # Other RES
    elif pemmdb_tech == "Other RES":
        df = pd.read_excel(
            fn,
            sheet_name="Other RES",
        )

    # Reserves
    elif pemmdb_tech == "Reserves":
        df = pd.read_excel(
            fn,
            sheet_name="Reserves",
        )

    # DSR
    elif pemmdb_tech == "DSR":
        df = pd.read_excel(
            fn,
            sheet_name="DSR",
        )

    # Battery
    elif pemmdb_tech == "Battery":
        df = pd.read_excel(
            fn,
            sheet_name="Battery",
        )

    # Electrolyser
    elif pemmdb_tech == "Electrolyser":
        df = pd.read_excel(
            fn,
            sheet_name="Electrolyser",
        )

    return df


def read_pemmdb_must_runs(
    node: str,
    pemmdb_dir: str,
    cyear: str,
    pyear: str,
    pemmdb_tech: str,
    sns: pd.DatetimeIndex,
) -> pd.Series:
    fn = Path(
        pemmdb_dir,
        str(pyear),
        f"PEMMDB_{node.replace('GB', 'UK')}_NationalTrends_{pyear}.xlsx",
    )

    if not os.path.isfile(fn):
        return None

    # Conventionals & Hydrogen
    if pemmdb_tech in CONVENTIONALS or pemmdb_tech == "Hydrogen":
        # read must runs
        must_runs = (
            pd.read_excel(
                fn,
                sheet_name="Thermal",
                skiprows=7,  # first seven rows are metadata
                usecols=[0, 1, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18],
                names=[
                    "carrier",
                    "type",
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
            .assign(
                **{
                    "carrier": lambda df: df.carrier.ffill(),
                    "type": lambda df: df.type.ffill().fillna(df.carrier),
                }
            )
            .query("unit == '% of installed capacity'")
            .drop(columns=["unit"])
            .set_index(["carrier", "type"])
            .div(1e2)  # percentage conversion
            .query("carrier == @pemmdb_tech")
        )
        must_runs = must_runs.set_index(must_runs.index.map("_".join)).T

        index_year = pd.date_range(
            start=f"{cyear}-01-01",
            periods=8760,  # 53 weeks
            freq="h",
        )
        df = pd.DataFrame({"month": index_year.month_name().str[:3]}, index=index_year)

        for type in must_runs.columns:
            df[type] = df["month"].map(must_runs[type])
        df.drop(columns="month", inplace=True)

        profile = (
            df.reset_index(names="time")
            .melt(id_vars="time", var_name="type", value_name="p_min_pu")
            .assign(carrier=pemmdb_tech, bus=node)
            .set_index(["time", "bus", "carrier", "type"])
            .to_xarray()
            .sel(time=sns)  # filter for snapshots only
        )

    # Wind
    elif pemmdb_tech == "Wind":
        df = pd.read_excel(fn, sheet_name="Wind")

    # Solar
    elif pemmdb_tech == "Solar":
        df = pd.read_excel(
            fn,
            sheet_name="Solar",
            # skiprows=1,
            # usecols=lambda name: name == "Day"
            # or name == "Week"
            # or name == "ShortName"
            # or name == "Variable"
            # or name == int(cyear),
            # sheet_name=f"{hydro_tech} - Year Dependent",
        )

    # Hydro
    elif pemmdb_tech == "Hydro":
        df = pd.read_excel(
            fn,
            sheet_name="Hydro",
            # skiprows=1,
            # usecols=lambda name: name == "Day"
            # or name == "Week"
            # or name == "ShortName"
            # or name == "Variable"
            # or name == int(cyear),
            # sheet_name=f"{hydro_tech} - Year Dependent",
        )

    # Other RES
    elif pemmdb_tech == "Other RES":
        df = pd.read_excel(
            fn,
            sheet_name="Other RES",
        )

    # Reserves
    elif pemmdb_tech == "Reserves":
        df = pd.read_excel(
            fn,
            sheet_name="Reserves",
        )

    # DSR
    elif pemmdb_tech == "DSR":
        df = pd.read_excel(
            fn,
            sheet_name="DSR",
        )

    # Battery
    elif pemmdb_tech == "Battery":
        df = pd.read_excel(
            fn,
            sheet_name="Battery",
        )

    # Electrolyser
    elif pemmdb_tech == "Electrolyser":
        df = pd.read_excel(
            fn,
            sheet_name="Electrolyser",
        )

    return profile


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_pemmdb_data",
            clusters="all",
            planning_horizons=2030,
            tech="Gas",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Climate year from snapshots
    sns = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)
    cyear = sns[0].year
    if int(cyear) < 1982 or int(cyear) > 2019:
        logger.warning(
            f"Snapshot year {cyear} doesn't match available TYNDP data. Falling back to 2009."
        )
        cyear = 2009

    # Planning year
    pyear = safe_pyear(
        snakemake.wildcards.planning_horizons,
        available_years=snakemake.params.available_years,
        source="PEMMDB",
    )

    # Parameters
    onshore_buses = pd.read_csv(snakemake.input.busmap, index_col=0)
    nodes = onshore_buses.index
    pemmdb_dir = snakemake.input.pemmdb_dir
    tech = str(snakemake.wildcards.tech)

    # Load and prep inflow data
    tqdm_kwargs_1 = {
        "ascii": False,
        "unit": " nodes",
        "total": len(nodes),
        "desc": "Loading PEMMDB capacities",
    }

    tqdm_kwargs_2 = {
        "ascii": False,
        "unit": " nodes",
        "total": len(nodes),
        "desc": "Loading PEMMDB must-runs",
    }

    func_caps = partial(
        read_pemmdb_capacities,
        pemmdb_dir=pemmdb_dir,
        cyear=cyear,
        pyear=pyear,
        pemmdb_tech=tech,
        sns=sns,
    )

    func_must_runs = partial(
        read_pemmdb_must_runs,
        pemmdb_dir=pemmdb_dir,
        cyear=cyear,
        pyear=pyear,
        pemmdb_tech=tech,
        sns=sns,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        pemmdb_capacities = list(tqdm(pool.imap(func_caps, nodes), **tqdm_kwargs_1))
        profiles = [
            profile
            for profile in tqdm(pool.imap(func_must_runs, nodes), **tqdm_kwargs_2)
            if profile is not None
        ]

    pemmdb_capacities_df = (
        pd.concat(pemmdb_capacities, axis=0)[
            ["bus", "carrier", "type", "p_nom", "country"]
        ]  # reorder columns
    )

    # merge p_min_pu profiles into one xarray dataset
    ds = xr.merge(profiles)

    pemmdb_capacities_df.to_csv(snakemake.output.pemmdb_capacities)
    ds.to_netcdf(snakemake.output.pemmdb_p_min_pu)
