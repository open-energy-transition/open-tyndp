# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds per-node TYNDP wheeling charges between the e-market and prosumer nodes.

TYNDP applies a wheeling charge between the e-market and prosumer nodes to
represent distribution grid costs; it is directional (only charged for flow
from the e-market to the prosumer node, not the other way around).

Inputs
------

- `data/tyndp/.../2026/Line-data/WHEELING_CHARGES.xlsx`: TYNDP 2026 wheeling
  charges, `Prosumer` sheet, one row per node.

Outputs
-------

- `resources/wheeling_charges_tyndp.csv`: per-node wheeling charges in
  EUR/MWh, indexed by node, with columns `e_market_to_prosumer` and
  `prosumer_to_e_market`.
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

COLUMN_MAP = {
    "NODE": "node",
    "WHEELING CHARGE eMARKET - PROSUMER [€/MWh]": "e_market_to_prosumer",
    "WHEELING CHARGE PROSUMER - eMARKET [€/MWh]": "prosumer_to_e_market",
}


def load_wheeling_charges(fn: str) -> pd.DataFrame:
    charges = pd.read_excel(fn, sheet_name="Prosumer", engine="calamine")
    charges = charges.rename(columns=COLUMN_MAP).set_index("node")
    charges.index = charges.index.str.replace("UK", "GB")
    return charges[["e_market_to_prosumer", "prosumer_to_e_market"]]


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_wheeling_charges",
            configfiles="config/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    wheeling_charges = load_wheeling_charges(snakemake.input.wheeling_charges)
    wheeling_charges.to_csv(snakemake.output.wheeling_charges, index=True)
