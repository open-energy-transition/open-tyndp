# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This script is used to clean TYNDP Scenario Building input data to be used in the PyPSA-Eur workflow. The `snapshot` year is used as climatic year (`cyear`). For DE and GA, it must be one of the following years: 1995, 2008 or 2009. For NT, it must be between 1982 and 2019. If the `snapshot` is not one of these years, then the demand is set to 2009 electricity demand (2009 being considered as the most representative of the three years).

Depending on the scenario, different planning years (`pyear`) are available. DE and GA are defined for 2030, 2040 and 2050. NT scenario is only defined for 2030 and 2040.
"""

import logging
import pandas as pd
from _helpers import configure_logging, set_scenario_config, get_snapshots

logger = logging.getLogger(__name__)


def extract_grid_data(h2_grid_raw, direction="1"):
    border_ctrys = h2_grid_raw.Border.str.split("-", expand=True)
    border_i = {
        "1": [0, 1],
        "2": [1, 0],
    }
    h2_grid = (
        h2_grid_raw
        .assign(country0=border_ctrys[border_i[direction][0]], country1=border_ctrys[border_i[direction][1]])
        .assign(p_nom=h2_grid_raw[f"Summary Direction {direction}"].mul(1e3))  # Gw to MW
    )[["country0", "country1", "p_nom"]]

    def make_index(c):
        return "H2 pipeline " + c.country0 + " -> " + c.country1

    h2_grid.index = h2_grid.apply(make_index, axis=1)
    return h2_grid

def load_and_clean_h2_grid(fn, scenario="GA", pyear=2030):
    """
    Load and clean H2 reference grid and format data.
    Returns both the cleaned reference grid and the interzonal connections as dataframes.
    """

    if scenario in ["DE", "GA"]:
        interzonal_raw = pd.read_excel(fn, sheet_name="Hydrogen_Interzonal")
        interzonal_raw.columns = interzonal_raw.columns.str.title()

        if int(pyear) not in [2030, 2035, 2040, 2045, 2050]:
            logger.warning("Planning horizon doesn't match available TYNDP data. "
                           "Falling back to closest available year between 2030 and 2050.")
            pyear = min(max(2030, 5 * round(pyear/5)), 2050)
        scenario_dict = {
            "GA": "Global Ambition",
            "DE": "Distributed Energy",
        }
        scenario = scenario_dict[scenario]
        interzonal_filtered = (
            interzonal_raw
            .query("Scenario == @scenario and Year == @pyear ")
        )

        interzonal_p = extract_grid_data(interzonal_filtered, direction="1")
        interzonal_reversed = extract_grid_data(interzonal_filtered, direction="2")
        interzonal = pd.concat([interzonal_p, interzonal_reversed]).sort_index()

    elif scenario == "NT":
            logger.info(
                "No interzonal capacities for 'NT' scenario. Saving empty file for interzonal capacities.")
            interzonal = pd.DataFrame()
    else:
        raise ValueError("Unknown scenario requested. Please, choose from 'GA', 'DE' or 'NT'.")

    h2_grid_raw = pd.read_excel(fn)
    h2_grid_p = extract_grid_data(h2_grid_raw, direction="1")
    h2_grid_reversed = extract_grid_data(h2_grid_raw, direction="2")
    h2_grid = pd.concat([h2_grid_p, h2_grid_reversed]).sort_index()

    return h2_grid, interzonal


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("clean_tyndp_h2_reference_grid",
                                   opts="",
                                   clusters="100",
                                   ll="v1.0",
                                   sector_opts="",
                                   planning_horizons="2030",
                                   )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    cf = snakemake.params.scenario
    scenario, pyear = cf.get("scenario", "DE"), cf.get("year", 2030)
    cyear = get_snapshots(snakemake.params.snapshots)[0].year

    # Load and prep H2 reference grid and interzonal pipeline capacities
    h2_grid, interzonal = load_and_clean_h2_grid(snakemake.input.tyndp_reference_grid, scenario=scenario, pyear=pyear)

    # Save prepped H2 grid and interzonal
    h2_grid.to_csv(snakemake.output.h2_grid_prepped)
    interzonal.to_csv(snakemake.output.interzonal_prepped)
