# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: MIT
"""
Loads and cleans the TYNDP capacity trajectories for a given TYNDP scenario.

Outputs
-------
Cleaned csv file with all TYNDP trajectories (`p_nom_min`, `p_nom_max`) as well as build out limitations in long format.

- ``resources/tyndp_trajectories.csv`` in long format.
"""

import logging

import pandas as pd

from scripts._helpers import (
    SCENARIO_DICT,
    configure_logging,
    map_tyndp_carrier_names,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_trajectories",
            clusters="all",
            planning_horizons=2030,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    fn = snakemake.input.trajectories
    tyndp_scenario = snakemake.params.tyndp_scenario
    column_names = {
        "NODE": "bus",
        "SCENARIO": "scenario",
        "TECHNOLOGY": "investment_dataset_carrier",
        "YEAR": "pyear",
        "MIN CAPACITY [MW]": "p_nom_min",
        "MAX CAPACITY [MW]": "p_nom_max",
    }

    # TODO: How to add buildout information if even necessary
    # Trajectories other than Nuclear are only used for DE and GA scenarios and after 2030
    trajectories_id = [] if tyndp_scenario == "NT" else ["All"]
    df = (
        pd.read_excel(fn, sheet_name="GLOBAL")
        .rename(column_names, axis="columns")
        .replace(SCENARIO_DICT, regex=True)
        .query(
            "scenario == @tyndp_scenario or (scenario in @trajectories_id and pyear > 2030)"
        )
    )

    carrier_mapping_df = (
        pd.read_csv(snakemake.input.carrier_mapping)[
            ["investment_dataset_carrier", "open_tyndp_carrier", "open_tyndp_index"]
        ]
    ).dropna()

    df = map_tyndp_carrier_names(df, carrier_mapping_df, ["investment_dataset_carrier"])

    df.to_csv(snakemake.output.tyndp_trajectories, index=False)
