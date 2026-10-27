# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP Scenario Building exogenous gas demand for Open-TYNDP.

Reads the annual methane final energy demand per EU27 country from the
``Methane | Total Energy Demand`` rows of the TYNDP 2026 Supply Tool's
``All data`` sheet, available for 2030, 2035, 2040 and 2050. This excludes
methane for power generation, SMR, pyrolysis, district heating, CHP and hybrid
heating. Countries outside the EU27 have no exogenous gas demand.

Inputs
------

- ``data/tyndp/archive/2026/Supply-Tool/``: TYNDP 2026 Supply Tool workbooks.

Outputs
-------

- ``gas_demand_tyndp_{planning_horizons}.csv``: Annual exogenous gas demand (MWh)
  per EU27 country for the specified planning horizon.
"""

import logging
from pathlib import Path

import pandas as pd
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def load_gas_demand(fn: Path, planning_horizon: int) -> pd.Series:
    """Load the annual methane final energy demand (MWh) per EU27 country."""
    data = pd.read_excel(
        fn,
        sheet_name="All data",
        usecols="A:AE",
        index_col=[0, 1, 3],
        engine="calamine",
    )
    demand = data.loc[[("Methane", "Total Energy Demand", planning_horizon)]].iloc[0]
    # TWh to MWh
    return demand.drop("Unit").astype(float).mul(1e6).rename("p_nom")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_gas_demand",
            configfiles="config/test/config.tyndp.yaml",
            planning_horizons=2035,
            run="NT",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    planning_horizon = int(snakemake.wildcards.planning_horizons)
    fn = Path(snakemake.input.supply_tool, snakemake.params.supply_tool_file)

    logger.info(f"Processing gas demand from {fn.name}")
    demand = load_gas_demand(fn, planning_horizon)

    demand.to_csv(snakemake.output.gas_demand, index=True)
