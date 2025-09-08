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


def _read_thermal_capacities(fn: Path, node: str, pemmdb_tech: str) -> pd.DataFrame:
    """
    Read and clean thermal (conventionals & hydrogen) capacities.
    """
    # read capacities (p_nom)
    df = (
        pd.read_excel(
            fn,
            sheet_name="Thermal",
            skiprows=7,
            usecols=[0, 1, 2],
            names=["carrier", "type", "p_nom"],
        )
        .iloc[3:]
        .dropna(how="all")
        .assign(
            carrier=lambda x: x.carrier.ffill(),
            type=lambda x: x.type.fillna(x.carrier),
            p_nom=lambda x: pd.to_numeric(x.p_nom, errors="coerce"),
            efficiency=1.0,  # capacities are in MWel and no efficiencies are given
            bus=node,
            country=node[:2],
        )
        .query("carrier == @pemmdb_tech")
        .reset_index(drop=True)
    )

    if df.empty:
        logger.info(
            f"No PEMMDB data matches climate year {cyear} for '{pemmdb_tech}' at {node}."
        )
        return None

    return df


def _read_other_nonres_capacities(
    fn: Path, node: str, cyear: str, pemmdb_tech: str
) -> pd.DataFrame:
    """
    Read and clean `Other Non-RES` profiles.
    """
    # read raw data and filter for price band columns
    df = pd.read_excel(
        fn,
        sheet_name="Other Non-RES",
        skiprows=7,
        index_col=1,
        nrows=9,
    ).filter(like="Price Band")

    if df.empty:
        logger.info(
            f"No PEMMDB data available for '{pemmdb_tech}' and climate year {cyear} at node {node}."
        )
        return None

    column_names = [
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

    df = (
        df.set_axis(column_names)
        .T.assign(
            carrier="Other Non-RES",
            bus=node,
            country=node[:2],
            cyear_start=lambda x: pd.to_numeric(x.cyear_start, errors="coerce"),
            cyear_end=lambda x: pd.to_numeric(x.cyear_end, errors="coerce"),
            p_nom=lambda x: pd.to_numeric(x.p_nom, errors="coerce"),
            units_count=lambda x: pd.to_numeric(x.units_count, errors="coerce"),
            price=lambda x: pd.to_numeric(x.price, errors="coerce"),
            efficiency=lambda x: pd.to_numeric(x.efficiency, errors="coerce"),
            co2_factor=lambda x: pd.to_numeric(x.co2_factor, errors="coerce"),
        )
        .query("cyear_start <= @cyear and cyear_end >= @cyear")
        .reset_index(drop=True)
    )

    if df.empty:
        logger.info(
            f"No PEMMDB data matches climate year {cyear} for '{pemmdb_tech}' at {node}."
        )
        return None

    return df


def read_pemmdb_capacities(
    node: str,
    pemmdb_dir: str,
    cyear: str,
    pyear: int,
    pemmdb_tech: str,
) -> pd.DataFrame:
    """
    Read and clean capacities from PEMMDB for a given technology, planning and climate year.

    Parameters
    ----------
    node : str
        Node name to read data for.
    pemmdb_dir : str
        Path to directory containing PEMMDB data.
    cyear : int
        Climate year to read data for.
    pyear : int
        Planning year to read data for. Can be fallback year to available data.
    pemmdb_tech : str
        PEMMDB technology to read data for.

    Returns
    -------
    df : pd.DataFrame
        Ddataframe containing NT capacities (p_nom) in long format for the given pemmdb technology.
    """
    fn = Path(
        pemmdb_dir,
        str(pyear),
        f"PEMMDB_{node.replace('GB', 'UK')}_NationalTrends_{pyear}.xlsx",
    )

    if not os.path.isfile(fn):
        logger.info(f"No PEMMDB data available for {node} in {pyear}.")
        return None

    try:
        # Conventionals & Hydrogen
        if pemmdb_tech in CONVENTIONALS or pemmdb_tech == "Hydrogen":
            return _read_thermal_capacities(fn, node, pemmdb_tech)

        # Other Non-RES
        elif pemmdb_tech == "Other Non-RES":
            return _read_other_nonres_capacities(fn, node, cyear, pemmdb_tech)

        # Wind
        elif pemmdb_tech == "Wind":
            pass  # placeholder

        # Solar
        elif pemmdb_tech == "Solar":
            pass  # placeholder

        # Hydro
        elif pemmdb_tech == "Hydro":
            pass  # placeholder

        # Other RES
        elif pemmdb_tech == "Other RES":
            pass  # placeholder

        # Reserves
        elif pemmdb_tech == "Reserves":
            pass  # placeholder

        # DSR
        elif pemmdb_tech == "DSR":
            pass  # placeholder

        # Battery
        elif pemmdb_tech == "Battery":
            pass  # placeholder

        # Electrolyser
        elif pemmdb_tech == "Electrolyser":
            pass  # placeholder

        else:
            return None

    except Exception as e:
        raise Exception(
            f"Error reading capacities for {pemmdb_tech} at {node} for climate year {cyear} and planning year {pyear}: {e}"
        )


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
        pemmdb_capacities = [
            caps
            for caps in tqdm(pool.imap(func, nodes), **tqdm_kwargs)
            if caps is not None
        ]

    if not pemmdb_capacities:
        raise Exception(
            f"No PEMMDB capacities available for '{tech}' with climate year {cyear} and planning year {pyear}. "
            f"Please specify different technology, climate year or planning year."
        )

    pemmdb_capacities_df = (
        pd.concat(pemmdb_capacities, axis=0)[
            ["carrier", "bus", "type", "p_nom", "efficiency", "country"]
        ]  # # select and order relevant columns
    ).reset_index(drop=True)

    pemmdb_capacities_df.to_csv(snakemake.output.pemmdb_capacities)
