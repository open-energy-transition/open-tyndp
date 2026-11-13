# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP Scenario Building exogenous gas demand for Open-TYNDP.

Reads the annual methane final energy demand per EU27 country from the
``Methane | Total Energy Demand`` rows of the TYNDP 2026 Supply Tool's
``All data`` sheet, available for 2030, 2035, 2040 and 2050. This row is the
final energy demand only: the Supply Tool country sheets add methane for power
generation, heat production (district heating, CHP and hybrid heating) and
SMR/pyrolysis on top of it to obtain the total methane demand (e.g. ``AT!D7 =
SUM(D8:D12)``), and the Scenario Report labels it "Final Energy Demand from ETM
(excl. DH; HHPs)". Power generation, SMR, pyrolysis and hybrid heating are
modelled endogenously. District heating and CHP heat are not and are added to
the exogenous demand as ``Methane | Primary Energy Demand for District
Heating`` plus the larger of ``Methane | Primary Energy Demand for CHP (Heat)``
and ``Methane for Heat | Methane``, following the heat production formula of
the Supply Tool country sheets (e.g. ``AT!L60``). Countries outside the EU27
have no exogenous gas demand.

Inputs
------

- ``data/tyndp/archive/2026/Supply-Tool/``: TYNDP 2026 Supply Tool workbooks.

Outputs
-------

- ``gas_demand_tyndp_{horizon}.csv``: Annual exogenous gas demand (MWh)
  per EU27 country for the specified planning horizon.
"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def load_gas_demand(fn: Path, planning_horizon: int) -> pd.Series:
    """
    Load the annual exogenous methane demand (MWh) per EU27 country.

    Sum of the final energy demand and the methane for district heating and
    CHP heat. Of the two CHP heat rows, the larger one is used per country, as
    in the Supply Tool country sheets.

    Parameters
    ----------
    fn : Path
        Path to the TYNDP 2026 Supply Tool workbook.
    planning_horizon : int
        Planning horizon to select from the ``Year`` column.

    Returns
    -------
    pd.Series
        Annual exogenous methane demand in MWh, indexed by country.
    """
    data = (
        pd.read_excel(
            fn,
            sheet_name="All data",
            usecols="A:AE",
            index_col=[0, 1, 3],
            engine="calamine",
        )
        .xs(planning_horizon, level=2)
        .drop(columns="Unit")
        .astype(float)
    )

    def row(etm: str, parameter: str) -> pd.Series:
        return data.loc[[(etm, parameter)]].squeeze(axis=0)

    final = row("Methane", "Total Energy Demand")
    district_heating = row("Methane", "Primary Energy Demand for District Heating")
    chp_heat = pd.concat(
        [
            row("Methane", "Primary Energy Demand for CHP (Heat)"),
            row("Methane for Heat", "Methane"),
        ],
        axis=1,
    ).max(axis=1)
    logger.info(
        f"Methane demand {planning_horizon} (TWh): final {final.sum():.1f}, "
        f"district heating {district_heating.sum():.1f}, CHP heat {chp_heat.sum():.1f}"
    )

    # TWh to MWh
    return (final + district_heating + chp_heat).mul(1e6).rename("p_nom")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_gas_demand",
            configfiles="config/test/config.tyndp.yaml",
            horizon=2035,
            run="NT",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    planning_horizon = int(snakemake.wildcards.horizon)
    fn = Path(snakemake.input.supply_tool, snakemake.params.supply_tool_file)

    logger.info(f"Processing gas demand from {fn.name}")
    demand = load_gas_demand(fn, planning_horizon)

    demand.to_csv(snakemake.output.gas_demand, index=True)
