# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds per-node TYNDP wheeling charges between the e-market and prosumer nodes.

TYNDP applies a wheeling charge between the e-market and prosumer nodes to
represent distribution grid costs. The charge is applied in both directions,
but the TYNDP 2026 data only prices the e-market to prosumer direction; the
prosumer to e-market direction is zero for every node.

The prosumer nodes are taken from the sheet names of the prosumer demand
files rather than from the wheeling charges file, since the two do not agree:
CH00 has a prosumer node in the TYNDP market model but no entry in the
wheeling charges file (its implied charge is zero there, verified against the
``Marginal Cost`` columns of the NT+ time series dashboards), while UK00 and
UKNI have an entry but no prosumer node. Prosumer nodes without an entry are
given a zero charge and reported in a warning.

Inputs
------

- `data/tyndp/.../2026/Line-data/WHEELING_CHARGES.xlsx`: TYNDP 2026 wheeling
  charges, `Prosumer` sheet, one row per node.
- `data/tyndp/.../2026/Demand/{pyear}/ELECTRICITY_PROSUMER {pyear}.xlsx`: TYNDP
  2026 prosumer demand, one sheet per prosumer node.

Outputs
-------

- `resources/wheeling_charges_tyndp.csv`: per-node wheeling charges in
  EUR/MWh, indexed by node, with columns `e_market_to_prosumer` and
  `prosumer_to_e_market`.
"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, format_bz_names, set_scenario_config
from scripts.sb.build_tyndp_demand import DEMAND_TYPE_MAP

logger = logging.getLogger(__name__)

COLUMN_MAP = {
    "NODE": "node",
    "WHEELING CHARGE eMARKET - PROSUMER [€/MWh]": "e_market_to_prosumer",
    "WHEELING CHARGE PROSUMER - eMARKET [€/MWh]": "prosumer_to_e_market",
}


def get_prosumer_nodes(fn: str) -> pd.Index:
    """
    Collect the TYNDP prosumer nodes from the prosumer demand files.

    Parameters
    ----------
    fn : str
        Path to the base directory containing per-year demand data
        subdirectories.

    Returns
    -------
    pd.Index
        Prosumer node names, formatted to Open-TYNDP's naming convention.

    Raises
    ------
    FileNotFoundError
        If no prosumer demand file is found.
    """
    prefix = DEMAND_TYPE_MAP["electricity_prosumer"]
    matches = sorted(Path(fn).glob(f"*/{prefix}[ _][0-9][0-9][0-9][0-9].xlsx"))

    if not matches:
        raise FileNotFoundError(f"No prosumer demand file found in {fn}")

    nodes = {
        sheet.removesuffix("_corrected")
        for match in matches
        for sheet in pd.ExcelFile(match, engine="calamine").sheet_names
    }
    return pd.Index(format_bz_names(pd.Series(sorted(nodes))))


def load_wheeling_charges(fn: str, nodes: pd.Index) -> pd.DataFrame:
    """
    Load the TYNDP wheeling charges for the given prosumer nodes.

    Parameters
    ----------
    fn : str
        Path to the TYNDP wheeling charges file.
    nodes : pd.Index
        Prosumer nodes to report charges for.

    Returns
    -------
    pd.DataFrame
        Wheeling charges in EUR/MWh indexed by prosumer node, with columns
        `e_market_to_prosumer` and `prosumer_to_e_market`.
    """
    charges = pd.read_excel(fn, sheet_name="Prosumer", engine="calamine")
    charges = charges.rename(columns=COLUMN_MAP).set_index("node")
    charges.index = format_bz_names(charges.index.to_series())
    charges = charges[["e_market_to_prosumer", "prosumer_to_e_market"]].reindex(nodes)

    missing = charges.index[charges.isna().any(axis=1)]
    if not missing.empty:
        logger.warning(
            f"No TYNDP wheeling charge given for {len(missing)} prosumer node(s), "
            f"assuming a zero charge for: {', '.join(missing)}"
        )

    return charges.fillna(0)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_wheeling_charges",
            configfiles="config/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    nodes = get_prosumer_nodes(snakemake.input.demand)
    wheeling_charges = load_wheeling_charges(snakemake.input.wheeling_charges, nodes)
    wheeling_charges.to_csv(snakemake.output.wheeling_charges, index=True)
