# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Applies the time aggregation, using snapshot weightings, on the sector-coupled network.

Description
-----------
Reads the snapshot weightings from the CSV file prepared in `build_snapshot_weightings`
and applies it on the time-varying network data prepared in `prepare_sector_network.py`.
Hourly ramp limits are scaled to the duration of each aggregated snapshot.
"""

import logging

import numpy as np
import pandas as pd
import pypsa

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)

logger = logging.getLogger(__name__)


def scale_ramp_limits(n: pypsa.Network) -> None:
    """
    Scale hourly ramp limits to match snapshot durations for generators and links.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network with ramp limits given per unit of nominal power per hour.

    Returns
    -------
    None
        Modifies the network in-place. Ramp limits that are constant across
        snapshots are set as static values. Otherwise, they are set as
        time-varying values and the static ones are reset to their default.
    """
    hours = n.snapshot_weightings.generators
    for c in n.components[{"Generator", "Link"}]:
        for attr in ["ramp_limit_up", "ramp_limit_down"]:
            limits = n.get_switchable_as_dense(c.name, attr).dropna(axis=1, how="all")
            if limits.empty:
                continue
            logger.info(
                f"Scale {attr} of {limits.shape[1]} {c.list_name} to match snapshot durations."
            )
            scaled = limits.mul(hours, axis=0).clip(upper=1.0)
            is_constant = scaled.eq(scaled.iloc[0]).all()
            c.static.loc[limits.columns, attr] = scaled.iloc[0].where(
                is_constant, c.defaults.at[attr, "default"]
            )
            c.dynamic[attr] = scaled.loc[:, ~is_constant]


def set_temporal_aggregation(
    n: pypsa.Network, resolution: str | bool, snapshot_weightings_fn: str
) -> pypsa.Network:
    """
    Aggregate time-varying data to the given snapshots.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network with hourly resolution.
    resolution : str | bool
        Temporal resolution specification.
    snapshot_weightings_fn : str
        Path to CSV file containing snapshot weightings for aggregation.

    Returns
    -------
    pypsa.Network
        Network with aggregated temporal resolution.
    """

    if not resolution:
        logger.info("No temporal aggregation. Using native resolution.")
        return n
    elif "sn" in resolution.lower():
        # Representative snapshots are dealt with directly
        sn = int(resolution[:-2])
        logger.info("Use every %s snapshot as representative", sn)
        n.set_snapshots(n.snapshots[::sn])
        n.snapshot_weightings *= sn
        scale_ramp_limits(n)
        return n
    else:
        # Otherwise, use the provided snapshots
        logger.info(
            f"Apply {resolution} temporal aggregation using snapshot weightings."
        )
        snapshot_weightings = pd.read_csv(
            snapshot_weightings_fn, index_col=0, parse_dates=True
        )

        # Define a series used for aggregation, mapping each hour in
        # n.snapshots to the closest previous timestep in
        # snapshot_weightings.index
        aggregation_map = (
            pd.Series(
                snapshot_weightings.index.get_indexer(n.snapshots), index=n.snapshots
            )
            .replace(-1, np.nan)
            .ffill()
            .astype(int)
            .map(lambda i: snapshot_weightings.index[i])
        )

        m = n.copy(snapshots=[])
        m.set_snapshots(snapshot_weightings.index)
        m.snapshot_weightings = snapshot_weightings

        # Aggregate all time-varying data.
        for c in n.components:
            pnl = getattr(m, c.list_name + "_t")
            for k, df in c.dynamic.items():
                if not df.empty:
                    if c.list_name == "stores" and k == "e_max_pu":
                        pnl[k] = df.groupby(aggregation_map).min()
                    elif c.list_name == "stores" and k == "e_min_pu":
                        pnl[k] = df.groupby(aggregation_map).max()
                    else:
                        pnl[k] = df.groupby(aggregation_map).mean()

        scale_ramp_limits(m)

        return m


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "temporal_aggregation",
            configfiles="config/test/config.tyndp.yaml",
            opts="",
            clusters="all",
            sector_opts="",
            planning_horizons="2030",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    n_h = pypsa.Network(snakemake.input.network)

    n = set_temporal_aggregation(
        n=n_h,
        resolution=snakemake.params.time_resolution,
        snapshot_weightings_fn=snakemake.input.snapshot_weightings,
    )

    n.export_to_netcdf(snakemake.output[0])
