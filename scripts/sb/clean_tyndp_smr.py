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
    SCENARIO_DICT,
    configure_logging,
    get_h2_zone_buses,
    safe_pyear,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def load_smr_data(
    fn: str, buses_h2_file: str, pyear: int, scenario: str
) -> pd.DataFrame:
    """
    Load and clean TYNDP SMR capacity, must run and CCS information.

    Parameters
    ----------
    fn : str
        Path to Excel file containing TYNDP SMR data.
    buses_h2_file : str
        Path to the TYNDP H2 buses CSV file, used to resolve each country's
        Z1 (or, absent a dedicated Z1 bus, Z2) hydrogen bus.
    pyear : int
        Planning horizon to read SMR data for.
    scenario : str
        TYNDP scenario to filter for.

    Returns
    -------
    pd.DataFrame
        Cleaned TYNDP SMR data with capacity, must run and CCS information.
    """

    column_dict = {
        "YEAR": "year",
        "SCENARIO": "scenario",
        "NODE": "bus",
        "CAPACITY [MW]": "p_nom",
        "HEAT RATE [GJ/MWh]": "heat_rate",
        "VO&M CHARGE [€/MWh]": "marginal_cost",
        "MUST-RUN UNITS": "must_run",
        "CCS": "ccs",
    }

    replace_dict = SCENARIO_DICT | {"UK": "GB"}

    # SMR is a Z1-role technology: use each country's dedicated Z1 bus, or
    # fall back to its Z2 bus if it has no separate Z1 zone (matches
    # `buses_h2_z1_effective` in `add_h2_topology_tyndp`)
    h2_zone_buses = get_h2_zone_buses(buses_h2_file)
    bus_z1_effective = h2_zone_buses.z1.combine_first(h2_zone_buses.z2)

    # Read data and rename
    smr = (
        pd.read_excel(fn)
        .rename(columns=column_dict)
        .replace(replace_dict)
        .query("year == @pyear and scenario == @scenario")
        .assign(
            bus=lambda df: df.bus.map(bus_z1_effective),
            carrier=lambda df: np.where(df.ccs, "SMR CC", "SMR"),
            p_min_pu=lambda df: np.where(df.must_run, 1, 0),
            efficiency=lambda df: 3.6 / df.heat_rate,  # convert to [MW_CH4/MW_H2]
            p_nom=lambda df: df.p_nom / df.efficiency,  # convert to [MW_CH4]
            unit="MW_CH4",
        )
        .drop(columns=["heat_rate", "marginal_cost", "must_run", "ccs", "efficiency"])
    )

    smr.index = smr.bus + " " + smr.carrier

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
    buses_h2_fn = snakemake.input.buses_h2
    scenario = snakemake.params.tyndp_scenario

    # Fallback for NT scenario
    if scenario == "NT":
        pyear = safe_pyear(pyear, [2030, 2040])

    # Load and prep SMR data
    smr = load_smr_data(
        fn=smr_fn, buses_h2_file=buses_h2_fn, pyear=pyear, scenario=scenario
    )

    # Save clean H2 SMR data
    smr.to_csv(snakemake.output.smr_prepped)
