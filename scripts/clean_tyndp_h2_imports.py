# SPDX-FileCopyrightText: : 2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This script is used to load and clean TYNDP H2 import data.
"""

import logging

import pandas as pd
from _helpers import (
    configure_logging,
    make_index,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def load_import_data(fn):
    """
    Load and clean H2 import potentials and marginal cost for pipeline and shipping
    Returns the cleaned data as dataframe.

    Parameters
    ----------
    fn : str
        Path to Excel file containing TYNDP H2 imports data.

    Returns
    -------
    pd.DataFrame
        The function returns cleaned TYNDP H2 import potentials and marginal cost.
    """

    imports = pd.read_excel(fn)

    rename_dict = {
        "YEAR": "Year",
        "SCENARIO": "Scenario",
        "CORRIDOR": "Corridor",
        "NODE FROM": "bus0",
        "NODE TO": "bus1",
        "MAX CAPACITY [MW]": "p_nom_link",
        "OFFER QUANTITY [MW]": "p_nom_generator",
        "OFFER PRICE [€/MWh]": "marginal_cost",
        "MAX ENERGY YEAR [GWh]": "e_sum_max",
    }

    scenario_dict = {
        "Distributed Energy": "DE",
        "Global Ambition": "GA",
        "National Trends": "NT",
    }

    imports = imports.rename(columns=rename_dict).replace(scenario_dict)
    imports.loc[:, "e_sum_max"] *= 1e3  # convert from GWh to MWh
    imports.loc[:, "Band"] = (
        imports.Corridor.str.split("-", expand=True).iloc[:, -1].str.lower()
    )
    imports.index = (
        imports.apply(make_index, axis=1, args=("H2 import",)) + " - " + imports.Band
    )

    return imports


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("clean_tyndp_h2_imports")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load and prep import potentials
    import_potentials = load_import_data(snakemake.input.import_potentials_raw)

    # Save prepped H2 import potentials
    import_potentials.to_csv(snakemake.output.import_potentials_prepped)
