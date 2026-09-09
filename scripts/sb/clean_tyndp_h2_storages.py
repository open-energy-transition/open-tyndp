# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script cleans and extracts the TYNDP H2 Storage data and saves it in a common pypsa friendly format.
"""

import logging

import numpy as np
import pandas as pd

from scripts._helpers import (
    SCENARIO_DICT,
    configure_logging,
    get_h2_zone_buses,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def load_h2_storage_data(
    fn: str, buses_h2_file: str, pyear: int, scenario: str
) -> pd.DataFrame:
    """
    Load and clean TYNDP H2 storage energy capacities as well as charge/discharge capacities and efficiencies.

    Parameters
    ----------
    fn : str
        Path to Excel file containing TYNDP H2 storage data.
    buses_h2_file : str
        Path to the TYNDP H2 buses CSV file, used to resolve each country's
        Z1/Z2 hydrogen bus.
    pyear : int
        Planning horizon to read H2 storage data for.
    scenario : str
        TYNDP scenario to filter for.

    Returns
    -------
    pd.DataFrame
        Cleaned TYNDP H2 storage data.
    """

    column_dict = {
        "YEAR": "year",
        "SCENARIO": "scenario",
        "NODE": "bus",
        "H2 ZONE": "h2_zone",
        "CAPACITY [GWh]": "e_nom",
        "MAX POWER [MW]": "p_nom_discharge",
        "MAX LOAD [MW]": "p_nom_charge",
        "MAX CAPACITY [GWh]": "e_nom_max",
        "MAX POWER EXPANSION [MW]": "p_nom_max_discharge",
        "MAX LOAD EXPANSION [MW]": "p_nom_max_charge",
        "CHARGE EFFICIENCY [%]": "efficiency_charge",
        "DISCHARGE EFFICIENCY [%]": "efficiency_discharge",
    }

    replace_dict = SCENARIO_DICT | {
        "UK": "GB",
        "ZONE 1": "H2 Z1",
        "ZONE 2": "H2 Z2",
        "All": "all",
    }

    # Tank storage is only added to countries with a dedicated Z1 bus (a
    # handful of countries); cavern storage is added to every country's Z2
    # bus (`add_h2_storage_tyndp`) - countries with no matching bus (e.g.
    # "ZONE 1" entries for countries without a Z1 zone, or "EU") are dropped
    h2_zone_buses = get_h2_zone_buses(buses_h2_file)

    # Read data and rename
    storages = (
        pd.read_excel(fn, sheet_name="TEMPLATE")
        .rename(columns=column_dict)
        .replace(replace_dict)
        .assign(
            e_nom=lambda df: df.e_nom * 1e3,  # [MWh]
            e_nom_max=lambda df: df.e_nom_max * 1e3,  # [MWh]
            efficiency_charge=lambda df: df.efficiency_charge / 100,  # [1]
            efficiency_discharge=lambda df: df.efficiency_discharge / 100,  # [1]
            storage_tech=lambda df: np.where(
                df.h2_zone == "H2 Z2", "cavern-storage", "tank-storage"
            ),
            bus=lambda df: np.where(
                df.h2_zone == "H2 Z2",
                df.bus.map(h2_zone_buses.z2),
                df.bus.map(h2_zone_buses.z1),
            ),
        )
        .dropna(subset=["bus"])
        .assign(bus=lambda df: df.bus + " " + df.storage_tech)
        .drop(columns="storage_tech")
    )

    # Manually fix 2030 expansion limits for NL
    # TODO: Remove if fixed
    err_entry_i = storages.query(
        "bus.str.contains('NL') and year == 2030 and h2_zone == 'H2 Z2'"
    ).index
    # scale up by missing decimal
    if (
        storages.at[err_entry_i.item(), "p_nom_max_charge"]
        < storages.at[err_entry_i.item(), "p_nom_charge"]
    ):
        storages.loc[err_entry_i, ["p_nom_max_charge"]] *= 10
    if (
        storages.at[err_entry_i.item(), "p_nom_max_discharge"]
        < storages.at[err_entry_i.item(), "p_nom_discharge"]
    ):
        storages.loc[err_entry_i, ["p_nom_max_discharge"]] *= 10

    storages = storages.loc[
        ((storages.scenario == scenario) | (storages.scenario == "all"))
        & (storages.year == pyear)
    ]

    return storages


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_tyndp_h2_storages",
            planning_horizons=2030,
            configfiles="config/test/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    pyear = int(snakemake.wildcards.planning_horizons)
    h2_storage_fn = snakemake.input.h2_storages
    buses_h2_fn = snakemake.input.buses_h2
    scenario = snakemake.params.tyndp_scenario

    # Load and prep H2 storage data
    h2_storages = load_h2_storage_data(
        fn=h2_storage_fn, buses_h2_file=buses_h2_fn, pyear=pyear, scenario=scenario
    )

    # Save clean H2 Storage data
    h2_storages.to_csv(snakemake.output.h2_storages_prepped, index=False)
