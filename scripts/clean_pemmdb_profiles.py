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
]


def read_pemmdb_profiles(
    node: str,
    pemmdb_dir: str,
    tyndp_scenario: str,
    cyear: str,
    pyear: int,
    pyear_i: int,
    pemmdb_tech: str,
    sns: pd.DatetimeIndex,
    index_year: pd.DatetimeIndex,
) -> xr.Dataset:
    fn = Path(
        pemmdb_dir,
        str(pyear),
        f"PEMMDB_{node.replace('GB', 'UK')}_NationalTrends_{pyear}.xlsx",
    )

    if not fn.is_file():
        logger.warning(f"No PEMMDB data available for {node} in {pyear}.")
        return None

    # Conventionals & Hydrogen
    if pemmdb_tech in CONVENTIONALS or pemmdb_tech == "Hydrogen":
        # read data
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

        df = pd.DataFrame({"month": index_year.month_name().str[:3]}, index=index_year)

        for type in must_runs.columns:
            df[type] = df["month"].map(must_runs[type])
        df.drop(columns="month", inplace=True)

        profiles = (
            df.reset_index(names="time")
            .melt(id_vars="time", var_name="type", value_name="p_min_pu")
            .assign(carrier=pemmdb_tech, bus=node, p_max_pu=1.0)
            .set_index(["time", "bus", "carrier", "type"])
            .to_xarray()
            .sel(time=sns)  # filter for snapshots only
        )

        # Remove must-runs for DE and GA after 2030 as to 2024 TYNDP Methodology report, p.37
        if tyndp_scenario != "NT" and pyear_i > 2030:
            profiles = profiles.assign(p_min_pu=0.0)

    # Other RES
    elif pemmdb_tech == "Other RES":
        must_runs = pd.read_excel(
            fn,
            sheet_name="Other RES",
        )

    # Other Non-RES
    elif pemmdb_tech == "Other Non-RES":
        # read data
        df = pd.read_excel(
            fn,
            sheet_name="Other Non-RES",
            skiprows=7,  # first seven rows are metadata
            index_col=1,
        ).rename(index={"Installed capacity \n(MW)": "p_nom", "PEMMDB type(s)": "type"})

        # create filter mask for column that matches the specified climate year
        # start_years = pd.to_numeric(df.iloc[7], errors='coerce')
        start_years = pd.to_numeric(df.loc["Start climate year", :], errors="coerce")
        # end_years = pd.to_numeric(df.iloc[8], errors='coerce')
        end_years = pd.to_numeric(df.loc["End climate year", :], errors="coerce")
        mask = (start_years <= cyear) & (cyear <= end_years)

        if not mask.any():
            logger.warning(
                f"No data available for '{pemmdb_tech}' and climate year {cyear} at node {node}."
            )
            return None

        df = (
            df.loc[:, mask].set_axis(  # filter for climate year
                ["p_max"], axis="columns"
            )  # rename the column
        )

        # extract capacity for per unit calculation and plant type
        cap = pd.to_numeric(df.loc["p_nom", :], errors="coerce").values[0]
        type = df.loc["type", :].values[0]

        profiles = (
            df.reset_index(drop=True)
            .loc[10:, :]
            .assign(
                p_nom=cap,
                p_max_pu=lambda df: df.p_max
                / df.p_nom,  # calculate per unit availability
                p_min_pu=0.0,
                time=index_year,
                bus=node,
                carrier=pemmdb_tech,
                type=type,
            )
            .set_index(["time", "bus", "carrier", "type"])[
                ["p_min_pu", "p_max_pu"]
            ]  # only keep per unit columns
            .to_xarray()
            .sel(time=sns)
        )

    else:
        return None

    return profiles


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
    if int(cyear) not in [1995, 2008, 2009]:
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

    # Load and prep data
    tqdm_kwargs = {
        "ascii": False,
        "unit": " nodes",
        "total": len(nodes),
        "desc": "Loading PEMMDB profiles",
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

    # merge pemmdb profiles into one xarray dataset and save
    ds = xr.merge(profiles)
    ds.to_netcdf(snakemake.output.pemmdb_profiles)
