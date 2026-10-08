# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Build the TYNDP 2026 synthetic fuel production links for a given planning horizon.

Reads the maximum H2 to synthetic fuel capacities (sheet ``H2 LIMITS``) and saves
one row per non-zero link from a national H2 Zone 2 bus to the EU27 ``e-liquids``
or ``sng`` bus, with its capacity in MW_H2.
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_synfuels", planning_horizons="2040", run="NT"
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    year = int(snakemake.wildcards.planning_horizons)

    links = (
        pd.read_excel(
            snakemake.input.synfuel_lines, sheet_name="H2 LIMITS", engine="calamine"
        )
        .query("YEAR == @year")
        .rename(
            columns={
                "NODE FROM": "bus0",
                "NODE TO": "bus1",
                "H2 - SYNTHETIC FUEL MAX CAPACITY [MW_h2]": "p_nom",
            }
        )
        .replace({"bus1": {"e_liquids": "e-liquids"}})
        .query("p_nom > 0")
        .loc[:, ["bus0", "bus1", "p_nom"]]
    )
    links.index = links.bus0 + "-" + links.bus1

    links.to_csv(snakemake.output.synfuel_links)
