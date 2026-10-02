# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds the TYNDP 2026 EV charging demand per fleet segment and the EV charging
station parameters.

TYNDP models EV charging per node in two locations, street charging on the
e-market node and home charging on the prosumer node, each with a fixed fleet
that charges on a fixed schedule and a flexible fleet. The raw EV demand from
`build_tyndp_demand.py` (demand types ``ev_market`` / ``ev_prosumer``, read from
``EV_FIXED_LOAD_PROFILES_ELECTRICITY_{MARKET,PROSUMER}``) covers the fixed fleet
only: divided by the charging efficiency it equals the ``{node} Street Fixed``
and ``{node} Prosumer Fixed`` charging loads reported in the NT+ time series
dashboards, whose element names are used for the output columns.

The raw EV demand is the energy stored in the vehicle batteries, i.e. before
charging losses. The charging efficiency is not applied here but exported
along with the other charging station parameters, to be applied on the
charging link.

Inputs
------

- `resources/demand_tyndp_ev_market_{planning_horizons}.csv`: Raw street EV
  charging demand per node, from `build_tyndp_demand.py`.
- `resources/demand_tyndp_ev_prosumer_{planning_horizons}.csv`: Raw home EV
  charging demand per node, from `build_tyndp_demand.py`.
- `data/tyndp/.../2026/Electric-Vehicle/EV_CHARGING STATIONS.xlsx`: EU-wide
  charging station parameters per segment and planning horizon.

Outputs
-------

- `resources/ev_demand_tyndp_{planning_horizons}.csv`: EV charging demand in
  MW per segment, with columns ``{node} {Street|Prosumer} Fixed``.
- `resources/ev_charging_stations_tyndp_{planning_horizons}.csv`: Charging
  station parameters indexed by segment (``Street``, ``Prosumer``), with
  charge/discharge rates in MW per station, charge/discharge efficiencies
  (p.u.), use of system charge (EUR/MWh) and stations per EV.
"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.sb.build_tyndp_demand import drop_zero_demand_columns

logger = logging.getLogger(__name__)

SEGMENTS = {"ev_market": "Street", "ev_prosumer": "Prosumer"}

CHARGING_STATIONS_COLUMNS = {
    "PROSUMER/STREET/FAST": "segment",
    "MAX CHARGE RATE [kW]": "charge_rate",
    "MAX DISCHARGE RATE [kW]": "discharge_rate",
    "CHARGE EFFICIENCY [%]": "charge_efficiency",
    "DISCHARGE EFFICIENCY [%]": "discharge_efficiency",
    "USE OF SYSTEM CHARGE [EUR/MWh]": "use_of_system_charge",
    "UNITS [-]": "units",
}


def get_charging_stations(fn: str, pyear: int) -> pd.DataFrame:
    """
    Get the EU-wide EV charging station parameters for a planning horizon.

    Parameters
    ----------
    fn : str
        Path to ``EV_CHARGING STATIONS.xlsx``.
    pyear : int
        Planning horizon.

    Returns
    -------
    pd.DataFrame
        Charging station parameters indexed by segment, with rates in MW per
        station and efficiencies in p.u.

    Raises
    ------
    ValueError
        If the file holds per-node parameters instead of EU-wide ones.
    """
    df = pd.read_excel(fn, engine="calamine")
    if (df["NODE"] != "All").any():
        raise ValueError(
            f"Per-node EV charging station parameters are not supported: {fn}"
        )

    stations = (
        df.loc[df["YEAR"] == pyear]
        .rename(columns=CHARGING_STATIONS_COLUMNS)
        .set_index("segment")[list(CHARGING_STATIONS_COLUMNS.values())[1:]]
    )
    if stations.empty:
        raise ValueError(f"No EV charging station parameters for {pyear} in {fn}")

    stations[["charge_rate", "discharge_rate"]] /= 1e3
    stations[["charge_efficiency", "discharge_efficiency"]] /= 100
    return stations


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_ev_demand",
            planning_horizons="2030",
            configfiles="config/test/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    pyear = int(snakemake.wildcards.planning_horizons)
    ev_dir = Path(snakemake.input.ev_modelling)

    segments = []
    for demand_type, segment in SEGMENTS.items():
        demand = pd.read_csv(
            snakemake.input[demand_type], index_col=0, parse_dates=True
        )
        segments.append(demand.add_suffix(f" {segment} Fixed"))
    ev_demand = drop_zero_demand_columns(pd.concat(segments, axis=1))
    charging_stations = get_charging_stations(
        ev_dir / "EV_CHARGING STATIONS.xlsx", pyear
    )

    ev_demand.to_csv(snakemake.output.ev_demand, index=True)
    charging_stations.to_csv(snakemake.output.charging_stations, index=True)
