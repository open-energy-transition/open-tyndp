# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script cleans the TYNDP 2024 Visualisation Platform data for benchmarking.

It reads and processes the data, applying correct indexes, headers, and naming conventions. The data
structure is converted to long format, and units are standardized to MW for power and MWh for energy.
"""

import logging

import country_converter as coco
import pandas as pd

from scripts._helpers import (
    SCENARIO_DICT,
    configure_logging,
    get_snapshots,
    set_scenario_config,
)
from scripts.clean_tyndp_benchmark import _convert_units

logger = logging.getLogger(__name__)

cc = coco.CountryConverter()
EU27_COUNTRIES = [
    country.replace("Czechia", "Czech Republic")
    for country in cc.EU27.name_short.to_list()
]
EU27_MAP = pd.Series(cc.EU27as("ISO2").ISO2.to_list(), index=EU27_COUNTRIES)

CARRIER_MAP = {
    "Gas CCGT": "methane",
    "Gas OCGT": "methane",
    "Gas conventional": "methane",
    "Hard coal": "coal + other fossil",
    "Hard coal biofuel": "coal + other fossil",
    "Heavy oil": "oil",
    "Hydrogen CCGT": "hydrogen",
    "Hydrogen FC": "hydrogen",
    "Light oil": "oil",
    "Lignite": "coal + other fossil",
    "Lignite biofuel": "biofuels",
    "Nuclear": "nuclear",
    "Oil shale": "oil",
    "Oil shale biofuel": "biofuels",
    "Other RES": "small scale res",
    "Others non-RES": "coal + other fossil",
    "Pondage": "hydro and pumped storage",
    "Reservoir": "hydro and pumped storage",
    "Run-of-River": "hydro and pumped storage",
    "Solar PV": "solar",
    "Solar Thermal": "solar thermal",
    "Wind Offshore": "wind offshore",
    "Wind Onshore": "wind onshore",
    "Large scale batteries": "battery",
    "Demand Side Response Explicit": "demand shedding",
    "Pump Storage - Closed Loop (turbine)": "hydro and pumped storage",
    "Pump Storage - Open Loop (turbine)": "hydro and pumped storage",
}


def get_elec_demand(
    elec_demand_fn: str,
    loss_factors_fn: str,
    scenario: str,
    cyear: int,
    unit_conversion: dict,
) -> pd.DataFrame:
    """
    Read and process the electricity demand data file. Transmission losses are deduced from the data to align with
    the Open-TYNDP scope.

    Parameters
    ----------
    elec_demand_fn : str
        Path to the electricity demand data file.
    loss_factors_fn : str
        Path to the loss factors file.
    scenario : str
        Name of the scenario being processed.
    cyear : int
        Climate year.
    unit_conversion : dict
        Dictionary of unit conversions.

    Returns
    -------
    pd.DataFrame
        Processed electricity demand data (MWh) with standardized format.
    """
    # Get loss factors
    loss_factors = (
        pd.read_csv(loss_factors_fn, index_col=0)
        .melt(var_name="year", ignore_index=False)
        .assign(year=lambda x: x["year"].astype("int"))
        .set_index("year", append=True)
    )

    # Read and process demand data
    df = (
        pd.read_excel(elec_demand_fn)
        .replace(SCENARIO_DICT, regex=True)
        .query(
            f"Scenario==@scenario and "
            f"Climate_Year=='CY{cyear}' and "
            f"Technology=='Native Demand (excl. pump load and battery charge)' and "
            f"Country in @EU27_COUNTRIES"
        )
        .assign(country_iso2=lambda x: x.Country.map(EU27_MAP))
        .rename(columns={"Year": "year", "Value": "value"})
        .set_index(["country_iso2", "year"])
    )
    source_unit = df.Unit_Name.iloc[0]
    df = df[["value"]]

    df /= 1 + loss_factors.reindex(df.index)

    df = _convert_units(df, source_unit, unit_conversion)

    data = (
        df.groupby("year")
        .value.sum()
        .reset_index()
        .assign(
            carrier="aggregated",
            scenario=f"TYNDP {scenario}",
            table="elec_demand",
        )
    )

    return data


def get_power_capacities(
    caps_fn: str, flex_fn: str, scenario: str, cyear: int, unit_conversion: dict
) -> pd.DataFrame:
    """
    Read and process the power capacity data files.

    Parameters
    ----------
    caps_fn : str
        Path to the power capacity data file.
    flex_fn : str
        Path to the flexible power capacity data file.
    scenario : str
        Name of the scenario being processed.
    cyear : int
        Climate year.
    unit_conversion : dict
        Dictionary of unit conversions.

    Returns
    -------
    pd.DataFrame
        Processed power capacity data with standardized format.
    """
    data = []
    for fn in [caps_fn, flex_fn]:
        df = (
            pd.read_excel(fn)
            .replace(SCENARIO_DICT, regex=True)
            .query(
                f"Scenario==@scenario and "
                f"Climate_Year=='CY{cyear}' and "
                f"Property_Name == 'Installed Capacity' and "
                f"Country in @EU27_COUNTRIES"
            )
            .rename(
                columns={"Year": "year", "Value": "value", "Category_Detail": "carrier"}
            )
            .assign(carrier=lambda x: x.carrier.map(CARRIER_MAP))
        )

        df = _convert_units(df, df.Unit_Name.iloc[0], unit_conversion)

        df = (
            df.groupby(["year", "carrier"])
            .value.sum()
            .reset_index()
            .assign(
                scenario=f"TYNDP {scenario}",
                table="power_capacity",
            )
        )
        data.append(df)

    data = (
        pd.concat(data)
        .groupby(["year", "carrier", "scenario", "table"])
        .sum()
        .reset_index()
    )

    return data


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("clean_tyndp_vp_data", run="NT")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    scenario = snakemake.params["scenario"]
    unit_conversion = snakemake.params["unit_conversion"]
    cyear = get_snapshots(snakemake.params.snapshots)[0].year

    logger.info("Reading Visualisation Platform data")

    elec_demand = get_elec_demand(
        elec_demand_fn=snakemake.input.elec_demand,
        loss_factors_fn=snakemake.input.loss_factors,
        scenario=scenario,
        cyear=cyear,
        unit_conversion=unit_conversion,
    )

    power_capacities = get_power_capacities(
        caps_fn=snakemake.input.elec_supplymix,
        flex_fn=snakemake.input.elec_flex,
        scenario=scenario,
        cyear=cyear,
        unit_conversion=unit_conversion,
    )

    # TODO Extend to include additional data from the Visualisation Platform

    df = pd.concat([elec_demand, power_capacities]).assign(source="TYNDP 2024 Vis Plfm")
    df.to_csv(snakemake.output[0], index=False)
