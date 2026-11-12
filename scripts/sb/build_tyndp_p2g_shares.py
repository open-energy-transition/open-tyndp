# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP split shares of electrolyser capacities between H2 Z2 zones.

PEMMDB gives one electrolyser capacity per electricity node, while some
countries have several H2 Z2 zones (such as FR, SK). The TYNDP P2G e-market file
gives the electrolyser capacity per electricity node and H2 zone pair.
PEMMDB remains the source of the capacities, with the split shares
coming from the P2G e-market file.

Inputs
------

- Input file from `.../Hydrogen/P2G_emarket_template_corrected_requirecheck.xlsx`:
  TYNDP 2026 P2G e-market capacities, one row per year and
  electricity node/H2 zone pair.

Outputs
-------

- `resources/p2g_shares_{horizon}.csv`: share of the
  electrolyser capacity of `bus0` that goes to `bus1`, with columns `bus0`,
  `bus1`, `p_nom` and `share`.
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.sb._helpers import format_bz_names

logger = logging.getLogger(__name__)

COLUMN_MAP = {
    "YEAR": "year",
    "NODE_FROM": "bus0",
    "NODE_TO": "bus1",
    "MAX_FLOW [MW]": "p_nom",
}


def load_p2g_shares(fn: str, planning_horizon: int) -> pd.DataFrame:
    """
    Load the TYNDP P2G e-market electrolyser capacities and calculate the
    split shares.

    Parameters
    ----------
    fn : str
        Path to the TYNDP P2G e-market file.
    planning_horizon : int
        Planning horizon to calculate the shares for.

    Returns
    -------
    pd.DataFrame
        Split shares with columns `bus0`, `bus1`, `p_nom` and `share`.
    """
    p2g = pd.read_excel(fn, sheet_name="TEMPLATE", engine="calamine")
    p2g = p2g.rename(columns=COLUMN_MAP).query("year == @planning_horizon")
    if p2g.empty:
        raise ValueError(f"No P2G e-market data for {planning_horizon}.")

    p2g[["bus0", "bus1"]] = format_bz_names(p2g[["bus0", "bus1"]])
    p2g["share"] = p2g.p_nom / p2g.groupby("bus0").p_nom.transform("sum")

    return p2g[["bus0", "bus1", "p_nom", "share"]]


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_p2g_shares",
            horizon=2030,
            configfiles="config/test/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    planning_horizon = int(snakemake.wildcards.horizon)
    p2g_shares = load_p2g_shares(snakemake.input.p2g_emarket, planning_horizon)
    p2g_shares.to_csv(snakemake.output.p2g_shares, index=False)
