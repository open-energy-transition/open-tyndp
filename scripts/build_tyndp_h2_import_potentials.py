# SPDX-FileCopyrightText: : 2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Filters TYNDP H2 import potentials and marginal cost for a specific TYNDP scenario and a given year wildcard.
The function saves a csv file with TYNDP H2 import potentials and marginal cost filtered for a specific TYNDP scenario
and for a given year.
"""

import logging

import pandas as pd
from _helpers import (
    configure_logging,
    set_scenario_config,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_h2_import_potentials",
            planning_horizons=2030,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    scenario = snakemake.params.tyndp_scenario
    year = int(snakemake.wildcards.planning_horizons)

    # Load prepped import potentials and filter
    fn = snakemake.input.import_potentials_prepped
    import_potentials = pd.read_csv(fn, index_col=0)
    import_potentials_filtered = import_potentials.query(
        "(Scenario == 'All' or Scenario == @scenario) and Year == @year"
    )

    # Save filtered H2 import potentials
    import_potentials_filtered.to_csv(snakemake.output.import_potentials_filtered)
