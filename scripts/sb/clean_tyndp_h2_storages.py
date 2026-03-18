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
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def load_h2_storage_data(
    fn: str, pyear: int, scenario: str, h2_zones_tyndp: bool
) -> pd.DataFrame:
    """
    Load and clean TYNDP H2 storage energy capacities as well as charge/discharge capacities and efficiencies.

    Parameters
    ----------
    fn : str
        Path to Excel file containing TYNDP H2 storage data.
    pyear : int
        Planning horizon to read H2 storage data for.
    scenario : str
        TYNDP scenario to filter for.
    h2_zones_tyndp : bool
        Whether H2 zonal split is modeled.

    Returns
    -------
    pd.DataFrame
        The function returns cleaned TYNDP H2 storage data.
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

    suffix = " H2 Z2" if h2_zones_tyndp else " H2"

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
            bus=lambda df: df.bus
            + np.where(df.h2_zone == "H2 Z2", suffix, " H2 Z1")
            + " "
            + np.where(df.h2_zone == "H2 Z2", "cavern-storage", "tank-storage"),
        )
    )

    # Manually fix 2030 expansion limits for NL
    # TODO: Remove if fixed
    err_entry_i = storages.query(
        "bus.str.contains('NL') and year == 2030 and h2_zone == 'H2 Z2'"
    ).index
    # scale up by missing decimal
    storages.loc[err_entry_i, ["p_nom_max_charge", "p_nom_max_discharge"]] *= 10

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
    scenario = snakemake.params.tyndp_scenario
    h2_zones_tyndp = snakemake.params.h2_zones_tyndp

    # Load and prep H2 storage data
    h2_storages = load_h2_storage_data(
        fn=h2_storage_fn, pyear=pyear, scenario=scenario, h2_zones_tyndp=h2_zones_tyndp
    )

    # Save clean H2 Storage data
    h2_storages.to_csv(snakemake.output.h2_storages_prepped, index=False)
