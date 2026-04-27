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
    convert_units,
    get_snapshots,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

cc = coco.CountryConverter()
EU27_COUNTRIES = [
    country.replace("Czechia", "Czech Republic")
    for country in cc.EU27.name_short.to_list()
]
EU27_MAP = pd.Series(cc.EU27as("ISO2").ISO2.to_list(), index=EU27_COUNTRIES)


# Mapping from TYNDP Visualization Platform (VP) carrier names to benchmark carrier names
def _load_vp_carrier_mapping(carrier_mapping_fn: str) -> dict:
    tech_map = pd.read_csv(carrier_mapping_fn)
    return (
        tech_map[
            [
                "tyndp_vp_carrier",
                "benchmarking_generation",
            ]
        ]
        .dropna(subset=["tyndp_vp_carrier"])
        .set_index("tyndp_vp_carrier")["benchmarking_generation"]
        .to_dict()
    )


def get_elec_demand(
    elec_demand_fn: str,
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

    df = convert_units(df, unit_col="Unit_Name")

    data = (
        df.groupby(["year", "unit"])
        .value.sum()
        .reset_index()
        .assign(
            carrier="final demand (inc. t&d losses, excl. pump storage )",
            scenario=f"TYNDP {scenario}",
            table="elec_demand",
            bus="EU27",
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
                columns={"Category_Detail": "carrier", "Year": "year", "Value": "value"}
            )
            .assign(carrier=lambda x: x.carrier.map(CARRIER_MAP))
        )

        df = convert_units(df, unit_col="Unit_Name")

        df = (
            df.groupby(["year", "carrier", "unit"])
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
        .groupby(["year", "carrier", "scenario", "table", "unit"])
        .sum()
        .reset_index()
        .assign(bus="EU27")
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

    # load carrier mapping
    CARRIER_MAP = _load_vp_carrier_mapping(snakemake.input.carrier_mapping)

    logger.info("Reading Visualisation Platform data")

    elec_demand = get_elec_demand(
        elec_demand_fn=snakemake.input.elec_demand,
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

    df = pd.concat([elec_demand, power_capacities]).assign(
        source="TYNDP 2024 Vis Pltfm"
    )
    df.to_csv(snakemake.output[0], index=False)
