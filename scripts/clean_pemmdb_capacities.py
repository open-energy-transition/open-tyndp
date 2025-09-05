# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT
"""
Loads and cleans the available PEMMDB v2.4 capacities from TYNDP data bundle for a given

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
Cleaned csv file with NT capacities (p_nom) in long format for the given pemmdb technology.
"""

import logging
import multiprocessing as mp
import os
from functools import partial
from pathlib import Path

import pandas as pd
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


def read_pemmdb_capacities(
    node: str,
    pemmdb_dir: str,
    cyear: str,
    pyear: int,
    pemmdb_tech: str,
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
                    "efficiency": 1.0,  # capacities are in MWel and no efficiencies are given
                    "type": lambda df: df["type"].fillna(df.carrier),
                    "bus": node,
                    "country": node[:2],
                }
            )
            .query("carrier == @pemmdb_tech")
        )

    # Other Non-RES
    elif pemmdb_tech == "Other Non-RES":
        df = (
            pd.read_excel(
                fn,
                sheet_name="Other Non-RES",
                skiprows=7,  # first seven rows are metadata
                index_col=1,
                nrows=9,  # capacity information is only at the top of the file
            )
            .filter(like="Price Band")
            .set_axis(
                [
                    "p_nom",
                    "units_count",
                    "type",
                    "purpose",
                    "price",
                    "efficiency",
                    "co2_factor",
                    "cyear_start",
                    "cyear_end",
                ]
            )
            .T.assign(
                **{
                    "carrier": "Other Non-RES",
                    "bus": node,
                    "country": node[:2],
                    "cyear_start": lambda x: x.cyear_start.astype(int),
                    "cyear_end": lambda x: x.cyear_end.astype(int),
                }
            )
            .query("cyear_start <= @cyear and cyear_end >= @cyear")
            .reset_index(drop=True)
        )

    # Wind
    elif pemmdb_tech == "Wind":
        df = pd.read_excel(
            fn,
            sheet_name="Wind",
        )

    # Solar
    elif pemmdb_tech == "Solar":
        df = pd.read_excel(
            fn,
            sheet_name="Solar",
        )

    # Hydro
    elif pemmdb_tech == "hydro":
        df = pd.read_excel(
            fn,
            sheet_name="Hydro",
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

    else:
        return None

    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_pemmdb_capacities",
            clusters="all",
            planning_horizons=2030,
            tech="Nuclear",
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
        "desc": "Loading PEMMDB capacities",
    }

    func = partial(
        read_pemmdb_capacities,
        pemmdb_dir=pemmdb_dir,
        cyear=cyear,
        pyear=pyear,
        pemmdb_tech=tech,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        pemmdb_capacities = list(tqdm(pool.imap(func, nodes), **tqdm_kwargs))

    pemmdb_capacities_df = (
        pd.concat(pemmdb_capacities, axis=0)[
            ["carrier", "bus", "type", "p_nom", "efficiency", "country"]
        ]  # # select and order relevant columns
    ).reset_index(drop=True)

    pemmdb_capacities_df.to_csv(snakemake.output.pemmdb_capacities)
