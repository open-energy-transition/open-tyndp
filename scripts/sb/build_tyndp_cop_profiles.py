# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Builds TYNDP 2026 heat pump coefficient of performance (COP) profiles.

TYNDP counterpart of `build_cop_profiles`: reads the hourly hybrid heat pump
COPs of the selected weather scenario and writes them in the upstream
``cop_profiles`` schema. Nodes with incomplete COPs or COPs outside
`COP_BOUNDS` are dropped.

Inputs
------

- `data/tyndp/.../2026/Heat/COP_{planning_horizons}.xlsx`

Outputs
-------

- `resources/cop_profiles_tyndp_{planning_horizons}.nc`
"""

import logging
from pathlib import Path

import pandas as pd
import xarray as xr

from scripts._helpers import (
    align_demand_to_snapshots,
    check_snapshot_year,
    configure_logging,
    get_snapshots,
    get_wscenario,
    set_scenario_config,
)
from scripts.sb.build_tyndp_demand import multiindex_to_datetimeindex

logger = logging.getLogger(__name__)

HEAT_SYSTEM = "hybrid"
HEAT_SOURCE = "air"
COP_BOUNDS = (1.0, 10.0)


def read_cop(
    fn: str, planning_horizon: int, wscenario: int, snapshots: pd.DatetimeIndex
) -> pd.DataFrame:
    """
    Read hourly COPs of all nodes for one planning horizon and weather scenario.

    Parameters
    ----------
    fn : str
        Path to the TYNDP 2026 ``Heat`` directory.
    planning_horizon : int
        Planning horizon selecting the ``COP_{planning_horizon}.xlsx`` file.
    wscenario : int
        Weather scenario column to read, e.g. 3 for ``WS003``.
    snapshots : pd.DatetimeIndex
        Snapshots to align the COPs to.

    Returns
    -------
    pd.DataFrame
        COPs indexed by ``time``, one column per node (``name``).
    """
    ws_code = f"WS{wscenario:03d}"
    cop_fn = Path(fn, f"COP_{planning_horizon}.xlsx")
    logger.info(f"Reading COPs for weather scenario {ws_code} from {cop_fn.name}")

    data = pd.read_excel(
        cop_fn,
        header=10,
        index_col=[0, 1],
        sheet_name=None,
        usecols=lambda name: name in ("Date", "Hour", ws_code),
        engine="calamine",
    )
    cop = pd.concat(data, axis=1).droplevel(1, axis=1)
    cop = multiindex_to_datetimeindex(cop, year=snapshots[0].year)
    cop.columns = cop.columns.str.replace("UK", "GB")
    cop = align_demand_to_snapshots(cop, snapshots)

    return cop.rename_axis(index="time", columns="name")


def drop_invalid_cop(cop: pd.DataFrame) -> pd.DataFrame:
    """
    Drop nodes with incomplete COPs or COPs outside `COP_BOUNDS`.

    Parameters
    ----------
    cop : pd.DataFrame
        COPs indexed by snapshot, one column per node.

    Returns
    -------
    pd.DataFrame
        `cop` restricted to nodes with complete and plausible COPs.
    """
    lower, upper = COP_BOUNDS
    invalid = cop.columns[cop.isna().any() | (cop.min() < lower) | (cop.max() > upper)]
    if not invalid.empty:
        logger.warning(
            f"Dropping nodes with incomplete or implausible COP data: {invalid.tolist()}"
        )

    return cop.drop(columns=invalid)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tyndp_cop_profiles",
            planning_horizons="2030",
            run="NT",
            configfiles="config/config.tyndp.yaml",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    planning_horizon = int(snakemake.wildcards.planning_horizons)
    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )
    check_snapshot_year(snapshots[0].year, snakemake.params.drop_leap_day)

    wscenario = get_wscenario(snakemake.params.wscenarios, planning_horizon)
    cop = read_cop(snakemake.input.heat, planning_horizon, wscenario, snapshots)
    cop = drop_invalid_cop(cop)

    xr.DataArray(cop).expand_dims(
        heat_system=[HEAT_SYSTEM], heat_source=[HEAT_SOURCE]
    ).to_netcdf(snakemake.output.cop_profiles)
