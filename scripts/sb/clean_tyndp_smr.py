# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script cleans and extracts the TYNDP SMR data and saves it in a common pypsa friendly format.
"""

import logging

import numpy as np
import pandas as pd

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def load_smr_data(fn: str, pyear: int, scenario: str) -> pd.DataFrame:
    """
    Load and clean TYNDP SMR capacity and CCS information.

    Parameters
    ----------
    fn : str
        Path to Excel file containing TYNDP SMR data.
    pyear : int
        Planning horizon to read SMR data for.
    scenario : str
        TYNDP scenario to filter for.

    Returns
    -------
    pd.DataFrame
        Cleaned TYNDP SMR data with capacity and CCS information.
    """

    column_dict = {
        "YEAR": "year",
        "SCENARIO": "scenario",
        "NODE": "bus",
        "CAPACITY [MW]": "p_nom",
        "HEAT RATE [GJ/MWh]": "heat_rate",
        "VO&M CHARGE [€/MWh]": "marginal_cost",
        "CCS": "ccs",
    }

    # TYNDP 2026 has no more scenario split: SCENARIO is always "All"
    replace_dict = {"All": "all"}

    # Read data and rename
    smr = (
        pd.read_excel(fn, sheet_name="TEMPLATE")
        .rename(columns=column_dict)
        .replace(replace_dict)
        .query("year == @pyear and (scenario == @scenario or scenario == 'all')")
        .assign(
            bus=lambda df: df.bus.str.replace("^UK", "GB", regex=True),
            carrier=lambda df: np.where(df.ccs, "SMR CC", "SMR"),
            # match the TYNDP market-output asset naming convention
            name_suffix=lambda df: np.where(df.ccs, "SMR CCS", "SMR"),
            p_min_pu=0,
            efficiency=lambda df: 3.6 / df.heat_rate,  # convert to [MW_CH4/MW_H2]
            p_nom=lambda df: df.p_nom / df.efficiency,  # convert to [MW_CH4]
            unit="MW_CH4",
        )
        .drop(columns=["heat_rate", "marginal_cost", "ccs", "efficiency"])
    )

    smr.index = smr.bus + " " + smr.name_suffix
    smr = smr.drop(columns="name_suffix")

    return smr


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "clean_tyndp_smr",
            planning_horizons=2030,
            configfiles="config/test/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    pyear = int(snakemake.wildcards.planning_horizons)
    smr_fn = snakemake.input.smr
    scenario = snakemake.params.tyndp_scenario

    # Load and prep SMR data
    smr = load_smr_data(fn=smr_fn, pyear=pyear, scenario=scenario)

    # Save clean H2 SMR data
    smr.to_csv(snakemake.output.smr_prepped)
