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
Cleaned csv file with NT capacities (p_nom) in long format for the given pemmdb technology. Also includes information on carrier, bus, type, efficiency, country and unit.
"""

import logging
import multiprocessing as mp
import os
from functools import partial
from pathlib import Path

import numpy as np
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

RENEWABLES = [
    "Solar",
    "Wind",
    "Hydro",
]

UNIT_CONVERSION = {
    "TWh": 1e6,  # TWh to MWh
    "GWh": 1e3,  # GWh to MWh
    "MWh": 1,  # MWh to MWh
    "kWh": 1e-3,  # kWh to MWh
    "GW": 1e3,  # GW to MW
    "MW": 1,  # MW to MW
    "kW": 1e-3,  # kW to MW
}


def _convert_units(
    df: pd.DataFrame,
    unit_conversion: dict[str, float],
    source_unit_col: str = "unit",
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
    unit_conversion : dict[str, float]
        Dictionary mapping units to conversion factors (to base unit).
    source_unit_col : str, default "unit
        Name of the column containing the source unit of the values.
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

    # Convert values using conversion factors
    conversion_factors = df[source_unit_col].map(unit_conversion)
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce") * conversion_factors

    # update unit column
    df["unit"] = df[source_unit_col].apply(
        lambda x: "MWh" if x in energy_units else "MW" if x in power_units else x
    )

    return df


def _read_thermal_capacities(
    fn: Path, node: str, cyear: int, pemmdb_tech: str
) -> pd.DataFrame:
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
            unit="MW",
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
    fn: Path, node: str, cyear: int, pemmdb_tech: str
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
            unit="MW",
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


def _parse_index_parts(
    index: pd.Index, pemmdb_tech: str
) -> tuple[pd.Series, pd.Series]:
    """
    Parse index parts of PEMMDB renewable capacity sheets based on technology type.
    """
    if pemmdb_tech == "Hydro":
        parts = index.str.split(r" - |capacity |\s*\(|\)", regex=True)
        type_series = parts.str[0]
        unit_series = parts.str[-2]
    else:
        parts = index.str.split(r"capacities |\s*\(|\)", regex=True)
        type_series = parts.str[1]
        unit_series = parts.str[2]

    return type_series, unit_series


def _read_res_capacities(
    fn: Path, node: str, cyear: int, pemmdb_tech: str, unit_conversion: dict[str, float]
) -> pd.DataFrame:
    """
    Read and clean `RES` (Solar, Wind, Hydro) capacities.
    """
    # read data
    df = pd.read_excel(
        fn, sheet_name=pemmdb_tech, skiprows=6, index_col=0, names=["p_nom"]
    ).dropna()

    if df.empty:
        logger.info(
            f"No PEMMDB data matches climate year {cyear} for '{pemmdb_tech}' at {node}."
        )
        return None

    # Induce type and unit from index
    types, units = _parse_index_parts(df.index, pemmdb_tech)

    df = df.assign(
        bus=node,
        country=node[:2],
        carrier=pemmdb_tech,
        efficiency=1.0,  # no efficiencies given but capacity in MWel
        type=types,
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
    df["type"] = np.where(
        df["element"] != "", df["type"] + " - " + df["element"], df["type"]
    )

    df = _convert_units(df, unit_conversion, "unit", "p_nom").reset_index(drop=True)

    return df


def _read_other_res_capacities(
    fn: Path, node: str, cyear: int, pemmdb_tech: str
) -> pd.DataFrame:
    """
    Read and clean `Other RES` capacities.
    """
    # read data
    df = (
        pd.read_excel(
            fn,
            sheet_name="Other RES",
            skiprows=6,
            usecols=[0, 1],
            index_col=0,
        )
        .iloc[[0]]
        .set_axis(["p_nom"], axis=1)
        .reset_index(drop=True)
        .assign(
            bus=node,
            country=node[:2],
            carrier=pemmdb_tech,
            efficiency=1.0,  # no efficiencies given but capacity in MWel
            type="Small Biomass, Geothermal, Marine, Waste and Not Defined",
            unit="MW",
        )
        .dropna()
    )

    if df.empty:
        logger.info(
            f"No PEMMDB data available for '{pemmdb_tech}' and climate year {cyear} at node {node}."
        )
        return None

    return df


def _read_electrolyser_capacities(
    fn: Path, node: str, cyear: int, pemmdb_tech: str
) -> pd.DataFrame:
    """
    Read and clean `Electrolyser` capacities.
    """
    # read data
    df = (
        pd.read_excel(
            fn,
            sheet_name=pemmdb_tech,
            skiprows=7,
            index_col=0,
        )
        .dropna(how="all", axis=0)
        .dropna(how="all", axis=1)
    )

    if df.empty:
        logger.info(
            f"No PEMMDB data available for '{pemmdb_tech}' and climate year {cyear} at node {node}."
        )
        return None

    column_names = [
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
        .assign(
            carrier=pemmdb_tech,
            bus=node,
            country=node[:2],
            type="Onshore grid connected",
            unit="MW",
        )
        .reset_index(drop=True)
    )

    return df


def _read_battery_capacities(
    fn: Path, node: str, cyear: int, pemmdb_tech: str
) -> pd.DataFrame:
    """
    Read and clean `Battery` capacities.
    """
    # read data
    df_raw = (
        pd.read_excel(fn, sheet_name="Battery", skiprows=7, index_col=0)
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
    types = ["Discharge", "Charge", "Store"]
    p_noms = pd.concat(
        [df_raw["p_nom_charge"], df_raw["p_nom_discharge"], df_raw["p_nom_store"]],
        axis=0,
    )

    df = pd.DataFrame(
        dict(
            p_nom=p_noms.values,
            efficiency=df_raw.efficiency[0],
            carrier=pemmdb_tech,
            bus=node,
            country=node[:2],
            type=types,
            unit=units,
        )
    )

    return df


def read_pemmdb_capacities(
    node: str,
    pemmdb_dir: str,
    cyear: int,
    pyear: int,
    pemmdb_tech: str,
    unit_conversion: dict[str, float],
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
    unit_conversion : dict[str, float]
        Dictionary mapping units to conversion factors (to base unit).

    Returns
    -------
    pd.DataFrame
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
            return _read_thermal_capacities(fn, node, cyear, pemmdb_tech)

        # Other Non-RES
        elif pemmdb_tech == "Other Non-RES":
            return _read_other_nonres_capacities(fn, node, cyear, pemmdb_tech)

        # Renewables (Solar, Wind, Hydro)
        elif pemmdb_tech in RENEWABLES:
            return _read_res_capacities(fn, node, cyear, pemmdb_tech, unit_conversion)

        # Other RES
        elif pemmdb_tech == "Other RES":
            return _read_other_res_capacities(fn, node, cyear, pemmdb_tech)

        # DSR
        elif pemmdb_tech == "DSR":
            pass  # placeholder

        # Battery
        elif pemmdb_tech == "Battery":
            return _read_battery_capacities(fn, node, cyear, pemmdb_tech)

        # Electrolyser
        elif pemmdb_tech == "Electrolyser":
            return _read_electrolyser_capacities(fn, node, cyear, pemmdb_tech)

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

    # Load and prep data
    tqdm_kwargs = {
        "ascii": False,
        "unit": " nodes",
        "total": len(nodes),
        "desc": f"Loading PEMMDB capacities for {tech}",
    }

    func = partial(
        read_pemmdb_capacities,
        pemmdb_dir=pemmdb_dir,
        cyear=cyear,
        pyear=pyear,
        pemmdb_tech=tech,
        unit_conversion=UNIT_CONVERSION,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        pemmdb_capacities = [
            caps
            for caps in tqdm(pool.imap(func, nodes), **tqdm_kwargs)
            if caps is not None
        ]

    if not pemmdb_capacities:
        logger.warning(
            f"No PEMMDB capacities available for '{tech}' with climate year {cyear} and planning year {pyear}. "
            f"Please specify different technology, climate year or planning year."
        )

    pemmdb_capacities_df = (
        pd.concat(pemmdb_capacities, axis=0)[
            ["carrier", "bus", "type", "p_nom", "efficiency", "country", "unit"]
        ]  # # select and order relevant columns
    ).reset_index(drop=True)

    pemmdb_capacities_df.to_csv(snakemake.output.pemmdb_capacities)
