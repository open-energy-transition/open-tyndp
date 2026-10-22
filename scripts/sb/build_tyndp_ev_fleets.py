# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds the TYNDP 2026 EV fleets per node: a fixed, a commuter and a
non-commuter fleet, each for street charging on the e-market node and home
charging on the prosumer node, named after the NT+ dashboard elements.

The fixed fleet demand is taken from `build_tyndp_demand.py`. The flexible
fleets are derived from the national EV stock and demand, with their driving,
availability and minimum state of charge profiles expanded on the TYNDP
reference calendar (365 days, starting on a Monday).

Inputs
------

- `resources/demand_tyndp_ev_{market,prosumer}_{planning_horizons}.csv`:
  fixed fleet demand per node.
- `data/tyndp/.../2026/Electric-Vehicle`: TYNDP 2026 EV modelling inputs.

Outputs
-------

- `resources/ev_demand_tyndp_{planning_horizons}.csv`: driving demand per
  fleet in MW, before charging losses.
- `resources/ev_fleets_tyndp_{planning_horizons}.csv`: fleet parameters, i.e.
  vehicles, grid side charging and V2G capacity, battery size, efficiencies
  and use of system charge.
- `resources/avail_profile_tyndp_{planning_horizons}.csv`: share of
  connected vehicles per flexible fleet.
- `resources/dsm_profile_tyndp_{planning_horizons}.csv`: minimum state of
  charge per flexible fleet.
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

SEGMENTS = {"ev_market": "Street", "ev_prosumer": "Prosumer"}

SEGMENT_SHARES = {
    "Street": "pEVs market node [%]",
    "Prosumer": "pEVs prosumer [%]",
}

V2G_SHARES = {"Street": "EV_V2G_Street", "Prosumer": "EV_V2G_Home"}

FLEXIBLE_FLEETS = {"Commuter": "COMMUTER", "Non-Commuter": "NON COMMUTER"}

CHARGING_STATIONS_COLUMNS = {
    "MAX CHARGE RATE [kW]": "charge_rate",
    "MAX DISCHARGE RATE [kW]": "discharge_rate",
    "CHARGE EFFICIENCY [%]": "charge_efficiency",
    "DISCHARGE EFFICIENCY [%]": "discharge_efficiency",
    "USE OF SYSTEM CHARGE [EUR/MWh]": "use_of_system_charge",
    "UNITS [-]": "units",
}

FLEET_COLUMNS = [
    "node",
    "segment",
    "fleet",
    "vehicles",
    "charge_p_nom",
    "discharge_p_nom",
    "e_nom",
    "charge_efficiency",
    "discharge_efficiency",
    "use_of_system_charge",
]


def read_parameters(
    fn: Path, pyear: int, year_column: str, index_column: str
) -> pd.DataFrame:
    """
    Read the EV parameters of a planning horizon from a TYNDP EV Excel file.

    Parameters
    ----------
    fn : Path
        Path to the Excel file.
    pyear : int
        Planning horizon.
    year_column : str
        Column holding the planning horizon.
    index_column : str
        Column to index by.

    Returns
    -------
    pd.DataFrame
        Parameters indexed by `index_column`.
    """
    df = pd.read_excel(fn, engine="calamine")
    df = df.loc[df[year_column] == pyear]
    return df.set_index(index_column, verify_integrity=True)


def get_charging_stations(fn: Path, pyear: int) -> pd.DataFrame:
    """
    Get the EU-wide EV charging station parameters for a planning horizon.

    Parameters
    ----------
    fn : Path
        Path to ``EV_CHARGING STATIONS.xlsx``.
    pyear : int
        Planning horizon.

    Returns
    -------
    pd.DataFrame
        Charging station parameters indexed by segment, with rates in MW per
        station and efficiencies in p.u.
    """
    stations = read_parameters(fn, pyear, "YEAR", "PROSUMER/STREET/FAST")
    stations = stations.rename(columns=CHARGING_STATIONS_COLUMNS)[
        list(CHARGING_STATIONS_COLUMNS.values())
    ]
    stations[["charge_rate", "discharge_rate"]] /= 1e3
    stations[["charge_efficiency", "discharge_efficiency"]] /= 100
    return stations


def expand_daily_profile(
    profiles: pd.DataFrame, column: str, index: pd.DatetimeIndex
) -> pd.Series:
    """
    Expand weekday and weekend profiles of a typical day on the snapshots,
    with the first day of the year being a Monday as in the TYNDP calendar.

    Parameters
    ----------
    profiles : pd.DataFrame
        Daily profiles indexed by hour from 1 to 24.
    column : str
        Column name with a ``{day_type}`` placeholder.
    index : pd.DatetimeIndex
        Snapshots.

    Returns
    -------
    pd.Series
        Profile values per snapshot.
    """
    weekday = profiles[column.format(day_type="WEEKDAY")].to_numpy()[index.hour]
    weekend = profiles[column.format(day_type="WEEKEND")].to_numpy()[index.hour]
    is_weekend = (index.dayofyear - 1) % 7 >= 5
    return pd.Series(np.where(is_weekend, weekend, weekday), index=index)


def get_fleets(fixed_demand: pd.DataFrame, ev_dir: Path, pyear: int) -> pd.DataFrame:
    """
    Get the parameters of the fixed and flexible EV fleets.

    Parameters
    ----------
    fixed_demand : pd.DataFrame
        Fixed fleet demand in MW.
    ev_dir : Path
        Path to the TYNDP 2026 EV modelling inputs.
    pyear : int
        Planning horizon.

    Returns
    -------
    pd.DataFrame
        Fleet parameters indexed by fleet, with the annual driving demand in
        MWh as ``energy``.
    """
    etm = read_parameters(ev_dir / "EV_ETM_DATA.xlsx", pyear, "YEAR", "COUNTRY")
    flex = read_parameters(ev_dir / "EV_FLEX_SHARE.xlsx", pyear, "YEAR", "Node")
    zone_shares = read_parameters(
        ev_dir / "EV_MARKET_ZONE_SHARES.xlsx", pyear, "year", "market_zone"
    )["split_factor"]
    misc = read_parameters(
        ev_dir / "EV_MISC_PARAMETERS.xlsx", pyear, "Year", "Year"
    ).loc[pyear]
    stations = get_charging_stations(ev_dir / "EV_CHARGING STATIONS.xlsx", pyear)

    nodes = fixed_demand.columns.str.split(" ").str[0].unique()
    base = pd.MultiIndex.from_product(
        [nodes, SEGMENTS.values()], names=["node", "segment"]
    ).to_frame(index=False)
    countries = base["node"].str[:2]
    share = base["segment"].map(SEGMENT_SHARES).map(misc) / 100
    share *= zone_shares.reindex(base["node"], fill_value=1).to_numpy()
    flex_share = flex.loc[base["node"], "EV_flex_share"].to_numpy() / 100
    base["v2g_share"] = [
        flex.at[node, V2G_SHARES[segment]] / 100
        for node, segment in zip(base["node"], base["segment"])
    ]
    base["vehicles"] = etm.loc[countries, "UNITS [#]"].to_numpy() * share
    base["energy"] = etm.loc[countries, "ELECTRITY DEMAND [TWh]"].to_numpy() * share
    base["energy"] *= 1e6
    base = base.join(stations, on="segment")

    fleet_shares = {"Fixed": 1 - flex_share} | dict.fromkeys(
        FLEXIBLE_FLEETS, flex_share / 2
    )
    fleets = pd.concat(
        [
            base.assign(
                fleet=fleet,
                vehicles=base["vehicles"] * fleet_share,
                energy=base["energy"] * fleet_share,
            )
            for fleet, fleet_share in fleet_shares.items()
        ],
        ignore_index=True,
    )
    fleets.index = fleets["node"] + " " + fleets["segment"] + " " + fleets["fleet"]
    fleets = fleets.loc[fleets["energy"] > 0]
    fixed = fleets["fleet"] == "Fixed"

    if len(fixed_demand) == 8760:
        deviation = fixed_demand.sum().div(fleets.loc[fixed, "energy"]) - 1
        mismatch = ~(deviation.abs() <= 0.01)
        if mismatch.any():
            raise ValueError(
                "Fixed EV fleet demand does not match the fixed share of the "
                f"national EV demand:\n{deviation[mismatch]}"
            )
    else:
        logger.info(
            "Skipping the fixed EV fleet demand reconciliation, as the demand "
            "does not cover the full year."
        )

    flexible_vehicles = fleets["vehicles"].where(~fixed, 0)
    fixed_p_nom = fixed_demand.max().reindex(fleets.index) / fleets["charge_efficiency"]
    fleets["charge_p_nom"] = (
        flexible_vehicles * fleets["units"] * fleets["charge_rate"]
    ).where(~fixed, fixed_p_nom)
    fleets["discharge_p_nom"] = (
        flexible_vehicles
        * fleets["v2g_share"]
        * fleets["units"]
        * fleets["discharge_rate"]
    )
    fleets["e_nom"] = flexible_vehicles * misc["Capacity [kWh]"] / 1e3

    return fleets


def get_flexible_profiles(
    flexible: pd.DataFrame, index: pd.DatetimeIndex, ev_dir: Path, pyear: int
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Get the driving demand, availability and minimum state of charge of the
    flexible EV fleets.

    Parameters
    ----------
    flexible : pd.DataFrame
        Flexible fleet parameters.
    index : pd.DatetimeIndex
        Snapshots.
    ev_dir : Path
        Path to the TYNDP 2026 EV modelling inputs.
    pyear : int
        Planning horizon.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        Driving demand in MW, availability and minimum state of charge in p.u.
    """
    driving_profiles = pd.read_excel(
        ev_dir / "DRIVING_PROFILES.xlsx", index_col="Hour", engine="calamine"
    )
    availability_profiles = pd.read_excel(
        ev_dir / "AVAILABILITY_PROFILES.xlsx", index_col="Hour", engine="calamine"
    )
    min_soc_profiles = pd.read_excel(
        ev_dir / "EV_MIN_SOC.xlsx",
        sheet_name="data",
        index_col="Hour",
        engine="calamine",
    )

    driving_column = "PASSENGER {day_type} [%]"
    year = pd.date_range(str(index.year[0]), periods=8760, freq="h")
    driving = expand_daily_profile(driving_profiles, driving_column, index)
    driving /= expand_daily_profile(driving_profiles, driving_column, year).sum()
    demand = pd.DataFrame(
        np.outer(driving, flexible["energy"]), index=index, columns=flexible.index
    )
    avail_profile = pd.DataFrame(
        {
            name: expand_daily_profile(
                availability_profiles,
                f"{segment.upper()} {FLEXIBLE_FLEETS[fleet]} {{day_type}} [%]",
                index,
            )
            for name, segment, fleet in flexible[["segment", "fleet"]].itertuples()
        }
    )
    dsm_profile = pd.DataFrame(
        {
            name: expand_daily_profile(
                min_soc_profiles, f"{segment.upper()} {{day_type}} {pyear} [%]", index
            )
            for name, segment in flexible["segment"].items()
        }
    )

    return demand, avail_profile, dsm_profile


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_ev_fleets",
            planning_horizons="2030",
            configfiles="config/test/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    pyear = int(snakemake.wildcards.planning_horizons)
    ev_dir = Path(snakemake.input.ev_modelling)

    fixed_demand = pd.concat(
        [
            pd.read_csv(
                snakemake.input[demand_type], index_col=0, parse_dates=True
            ).add_suffix(f" {segment} Fixed")
            for demand_type, segment in SEGMENTS.items()
        ],
        axis=1,
    )
    fleets = get_fleets(fixed_demand, ev_dir, pyear)
    fixed = fleets.index[fleets["fleet"] == "Fixed"]
    flexible = fleets.loc[fleets["fleet"] != "Fixed"]
    flexible_demand, avail_profile, dsm_profile = get_flexible_profiles(
        flexible, fixed_demand.index, ev_dir, pyear
    )
    ev_demand = pd.concat([fixed_demand[fixed], flexible_demand], axis=1)

    ev_demand.to_csv(snakemake.output.ev_demand, index=True)
    fleets[FLEET_COLUMNS].to_csv(snakemake.output.ev_fleets, index=True)
    avail_profile.to_csv(snakemake.output.avail_profile, index=True)
    dsm_profile.to_csv(snakemake.output.dsm_profile, index=True)
