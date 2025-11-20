# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
This script loads and cleans the TYNDP H2 reference grid and interzonal connections for a given wildcard planning horizon
and TYNDP scenario as defined in the config file. The reference grid contains data for the TYNDP planning year 2030,
while depending on the scenario, different planning years (`pyear`) are available for the interzonal connections.
DE and GA are defined for 2030, 2035, 2040, 2045 and 2050. For the NT scenario no interzonal capacities are defined.
"""

import logging

import numpy as np
import pandas as pd

from scripts._helpers import (
    SCENARIO_DICT,
    configure_logging,
    extract_grid_data_tyndp,
    get_snapshots,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def load_h2_interzonal_connections(fn, scenario="GA", pyear=2030):
    """
    Load and clean H2 interzonal connections.
    Returns the cleaned interzonal connections as dataframe.

    Parameters
    ----------
    fn : str
        Path to Excel file containing H2 interzonal data.
    scenario : str
        TYNDP scenario to use for interzonal connection data.
        Possible options are:
        - 'GA'
        - 'DE'
        - 'NT'
    pyear : int
        TYNDP planning horizon to use for interzonal connection data.
        Possible options are:
        - 2030
        - 2035
        - 2040
        - 2045
        - 2050

    Returns
    -------
    pd.DataFrame
        The function returns cleaned TYNDP H2 interzonal connections.
    """

    if scenario in ["DE", "GA"]:
        interzonal_raw = pd.read_excel(fn, sheet_name="Hydrogen_Interzonal")
        interzonal_raw.columns = interzonal_raw.columns.str.title()

        if int(pyear) not in [2030, 2035, 2040, 2045, 2050]:
            logger.warning(
                "Planning horizon doesn't match available TYNDP data. "
                "Falling back to closest available year between 2030 and 2050."
            )
            pyear = np.clip(5 * (pyear // 5), 2030, 2050)
        scenario_dict = {v: k for k, v in SCENARIO_DICT.items()}
        scenario = scenario_dict[scenario]
        interzonal_filtered = interzonal_raw.query(
            "Scenario == @scenario and Year == @pyear "
        )

        interzonal = extract_grid_data_tyndp(interzonal_filtered, "H2 pipeline")
        # convert from GW to PyPSA base unit MW as raw H2 reference grid data is given in GW
        interzonal["p_nom"] = interzonal.p_nom.mul(1e3)

    elif scenario == "NT":
        logger.info(
            "No interzonal capacities for 'NT' scenario. Saving empty file for interzonal capacities."
        )
        interzonal = pd.DataFrame()
    else:
        raise ValueError(
            "Unknown scenario requested. Please, choose from 'GA', 'DE' or 'NT'."
        )

    return interzonal


def load_h2_grid(fn):
    """
    Load and clean H2 reference grid and format data.
    Returns the cleaned reference grid as dataframe.

    Parameters
    ----------
    fn : str
        Path to Excel file containing H2 reference grid data.

    Returns
    -------
    pd.DataFrame
        The function returns the cleaned TYNDP H2 reference grid.
    """

    h2_grid_raw = pd.read_excel(fn)
    h2_grid = extract_grid_data_tyndp(h2_grid_raw, "H2 pipeline")
    # convert from GW to PyPSA base unit MW as raw H2 reference grid data is given in GW
    h2_grid["p_nom"] = h2_grid.p_nom.mul(1e3)

    # remove H2 pipeline capacity between CH and ITIB as of updated model topology (https://tyndp2024.entsog.eu/h2igi-report/)
    ch_ibit_i = h2_grid.query(
        "index.str.contains('CH') and index.str.contains('IT')"
    ).index
    h2_grid.loc[ch_ibit_i, "p_nom"] = 0.0

    return h2_grid


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_h2_network",
            planning_horizons=2030,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    scenario = snakemake.params.scenario
    pyear = int(snakemake.wildcards.planning_horizons)
    cyear = get_snapshots(snakemake.params.snapshots)[0].year

    # Load and prep H2 reference grid and interzonal pipeline capacities
    h2_grid = load_h2_grid(fn=snakemake.input.tyndp_reference_grid)
    interzonal = load_h2_interzonal_connections(
        fn=snakemake.input.tyndp_reference_grid, scenario=scenario, pyear=pyear
    )

    # Save prepped H2 grid and interzonal
    h2_grid.to_csv(snakemake.output.h2_grid_prepped)
    interzonal.to_csv(snakemake.output.interzonal_prepped)
