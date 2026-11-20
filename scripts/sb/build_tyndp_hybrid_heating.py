# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds the TYNDP 2026 hybrid heat pump boiler topology.

Reads the CH4 and H2 hybrid heat pump boilers of the planning horizon from
``HHP_{CH4,H2}_BOILERS.xlsx`` and maps them to Open-TYNDP bus names. Each row
is one boiler with its hybrid heat pump node, the electricity bus of the node
and, for H2, the H2 bus it draws from. A node appears several times if its
boiler draws from several H2 buses. ``SK0E`` and ``SK0W`` have no electricity
bus of their own and are connected to ``SK00``. CH4 boilers have no
``fuel_bus``, as they draw from the gas bus of the node.

Inputs
------

- ``data/tyndp/.../2026/Heat/HHP_{CH4,H2}_BOILERS.xlsx``

Outputs
-------

- ``resources/hybrid_heating_boilers_{horizon}.csv``
"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.sb._helpers import format_bz_names

logger = logging.getLogger(__name__)


def read_boilers(heat_dir: str, fuel: str, planning_horizon: int) -> pd.DataFrame:
    """
    Read the hybrid heat pump boilers of one backup fuel and planning horizon.

    Parameters
    ----------
    heat_dir : str
        Path to the TYNDP 2026 ``Heat`` directory.
    fuel : str
        Backup fuel of the hybrid heat pump, ``"CH4"`` or ``"H2"``.
    planning_horizon : int
        Planning horizon to select from the ``YEAR`` column.

    Returns
    -------
    pd.DataFrame
        Boilers with columns ``fuel``, ``node``, ``bus`` and ``fuel_bus``.
    """
    boilers = pd.read_excel(
        Path(heat_dir, f"HHP_{fuel}_BOILERS.xlsx"), engine="calamine"
    ).query("YEAR == @planning_horizon")
    node = format_bz_names(boilers.NODE)

    return pd.DataFrame(
        {
            "fuel": fuel,
            "node": node,
            "bus": node.replace({"SK0E": "SK00", "SK0W": "SK00"}),
            "fuel_bus": boilers.NODE_FROM if fuel == "H2" else None,
        }
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_hybrid_heating",
            horizon="2030",
            run="NT",
            configfiles="config/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    planning_horizon = int(snakemake.wildcards.horizon)
    boilers = pd.concat(
        [
            read_boilers(snakemake.input.heat, fuel, planning_horizon)
            for fuel in ("CH4", "H2")
        ]
    )

    boilers.to_csv(snakemake.output.boilers, index=False)
