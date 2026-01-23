# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Solve CBA reference network with perfect foresight to extract marginal storage values (MSV).

This script solves the reference network for the full year (or at lower temporal resolution)
with fixed capacities and perfect foresight. The resulting dual variables from the store
energy balance constraints (mu_energy_balance) represent the marginal storage values that
capture seasonal patterns.

These MSV values are then used by solve_cba_network.py to guide rolling horizon dispatch
for both the reference network and all project networks.

References
----------
- PyPSA documentation on storage optimization
"""

import logging

import numpy as np
import pandas as pd
import pypsa

from scripts._benchmark import memory_logger
from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.solve_network import (
    add_co2_atmosphere_constraint,
    prepare_network,
)

logger = logging.getLogger(__name__)

# Default threshold for warning about extreme MSV values (EUR/MWh)
DEFAULT_MSV_EXTREME_THRESHOLD = 1000.0


def resample_network_temporal(n: pypsa.Network, resolution: str) -> pypsa.Network:
    """
    Resample network to a lower temporal resolution, including all time-series data.

    Follows the pattern from prepare_network.average_every_nhours() and
    prepare_sector_network.set_temporal_aggregation().

    Parameters
    ----------
    n : pypsa.Network
        Network to resample
    resolution : str
        Target resolution (e.g., "3H", "6H", "24H")

    Returns
    -------
    pypsa.Network
        New network with resampled snapshots and time-varying data
    """
    offset = resolution.lower()
    logger.info(f"Resampling network to {offset} resolution")

    # Create a new network with empty snapshots
    m = n.copy(snapshots=[])

    # Resample years separately to handle non-contiguous years
    # (pattern from time_aggregation.py)
    years = pd.DatetimeIndex(n.snapshots).year.unique()
    snapshot_weightings = []

    for year in years:
        sws_year = n.snapshot_weightings[n.snapshots.year == year]
        sws_year = sws_year.resample(offset).sum()
        snapshot_weightings.append(sws_year)

    snapshot_weightings = pd.concat(snapshot_weightings)

    # Drop rows with zero weight (time steps not in original snapshots)
    zeros_i = snapshot_weightings.query("objective == 0").index
    snapshot_weightings.drop(zeros_i, inplace=True)

    m.set_snapshots(snapshot_weightings.index)
    m.snapshot_weightings = snapshot_weightings

    # Resample ALL time-varying data
    # (pattern from prepare_sector_network.set_temporal_aggregation)
    for c in n.iterate_components():
        pnl = getattr(m, c.list_name + "_t")
        for k, df in c.pnl.items():
            if not df.empty:
                # Special handling for store energy limits
                if c.list_name == "stores" and k == "e_max_pu":
                    pnl[k] = df.resample(offset).min()
                elif c.list_name == "stores" and k == "e_min_pu":
                    pnl[k] = df.resample(offset).max()
                else:
                    # Default: use mean for averaging time-series
                    pnl[k] = df.resample(offset).mean()

    logger.info(
        f"Resampled from {len(n.snapshots)} to {len(m.snapshots)} snapshots "
        f"at {resolution} resolution"
    )

    return m


def validate_msv(
    n: pypsa.Network,
    seasonal_carriers: list[str] | None = None,
    extreme_threshold: float = DEFAULT_MSV_EXTREME_THRESHOLD,
) -> dict[str, dict]:
    """
    Validate MSV values and log warnings for extreme values.

    Parameters
    ----------
    n : pypsa.Network
        Network with solved MSV (mu_energy_balance)
    seasonal_carriers : list[str], optional
        List of carriers to validate. If None, validates all carriers.
    extreme_threshold : float, optional
        Threshold for warning about extreme MSV values (EUR/MWh).
        Default is 1000.

    Returns
    -------
    dict[str, dict]
        Dictionary with MSV statistics per carrier
    """
    stats = {}

    if n.stores_t.mu_energy_balance.empty:
        logger.warning("No mu_energy_balance found in stores")
        return stats

    carriers = seasonal_carriers or n.stores.carrier.unique().tolist()

    for carrier in carriers:
        stores_i = n.stores[n.stores.carrier == carrier].index
        stores_i = stores_i[stores_i.isin(n.stores_t.mu_energy_balance.columns)]

        if stores_i.empty:
            continue

        msv = n.stores_t.mu_energy_balance[stores_i]

        carrier_stats = {
            "count": len(stores_i),
            "min": msv.min().min(),
            "max": msv.max().max(),
            "mean": msv.mean().mean(),
            "std": msv.std().mean(),
        }
        stats[carrier] = carrier_stats

        # Log statistics
        logger.info(
            f"MSV for '{carrier}' ({len(stores_i)} stores): "
            f"min={carrier_stats['min']:.2f}, "
            f"max={carrier_stats['max']:.2f}, "
            f"mean={carrier_stats['mean']:.2f}, "
            f"std={carrier_stats['std']:.2f} EUR/MWh"
        )

        # Warn about extreme values
        if (msv.abs() > extreme_threshold).any().any():
            extreme_count = (msv.abs() > extreme_threshold).sum().sum()
            logger.warning(
                f"Extreme MSV values (|MSV| > {extreme_threshold} EUR/MWh) "
                f"detected for '{carrier}': {extreme_count} occurrences. "
                "This may indicate model issues or unusual scarcity patterns."
            )

    return stats


def extra_functionality_msv(
    n: pypsa.Network,
    snapshots: pd.DatetimeIndex,
) -> None:
    """
    Add constraints for MSV extraction solve.

    Applies a subset of constraints from the main solve to ensure consistency
    between MSV extraction and rolling horizon dispatch. Currently includes:
    - CO2 atmosphere constraint

    Note: Operational reserve margin is excluded as it's more relevant for
    short-term dispatch than for extracting seasonal storage values.

    Parameters
    ----------
    n : pypsa.Network
        Network with config attribute attached
    snapshots : pd.DatetimeIndex
        Simulation timesteps
    """
    add_co2_atmosphere_constraint(n, snapshots)


def solve_for_msv(
    n: pypsa.Network,
    solver_name: str = "highs",
    solver_options: dict | None = None,
    io_api: str | None = None,
    extra_functionality=None,
) -> tuple[str, str]:
    """
    Solve network with perfect foresight to extract MSV.

    Parameters
    ----------
    n : pypsa.Network
        Network to solve (with fixed capacities)
    solver_name : str
        Solver to use
    solver_options : dict, optional
        Solver options
    io_api : str, optional
        I/O API to use (e.g., "direct" for better performance)
    extra_functionality : callable, optional
        Function to add extra constraints (e.g., CO2 atmosphere constraint)

    Returns
    -------
    tuple[str, str]
        Status and termination condition
    """
    if solver_options is None:
        solver_options = {}

    logger.info(f"Solving network with {len(n.snapshots)} snapshots for MSV extraction")

    optimize_kwargs = {
        "solver_name": solver_name,
        "solver_options": solver_options,
        "assign_all_duals": True,  # Critical: needed to get mu_energy_balance
    }

    # Add io_api if specified (can provide ~20% speedup with "direct")
    if io_api:
        optimize_kwargs["io_api"] = io_api
        logger.info(f"Using io_api={io_api}")

    # Add extra_functionality for constraints (e.g., CO2 atmosphere)
    if extra_functionality:
        optimize_kwargs["extra_functionality"] = extra_functionality

    status, condition = n.optimize(**optimize_kwargs)

    if status == "ok":
        logger.info(f"Optimization successful: {condition}")
    else:
        logger.warning(f"Optimization failed: {status}, {condition}")

    return status, condition


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_cba_msv_extraction",
            planning_horizons="2030",
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    # Load reference network
    n = pypsa.Network(snakemake.input.network)
    logger.info(f"Loaded network with {len(n.snapshots)} snapshots")

    # Get CBA solving parameters (includes load_shedding: true)
    solving = snakemake.params.cba_solving
    solver_name = solving["solver"]["name"]
    solver_options_key = solving["solver"].get("options", None)
    solver_options = (
        solving["solver_options"].get(solver_options_key, {})
        if solver_options_key
        else {}
    )

    # Get io_api setting for performance optimization
    io_api = solving.get("options", {}).get("io_api", None)

    # Get MSV extraction settings
    msv_resolution = snakemake.params.msv_resolution

    # Resample to lower resolution if specified
    # This now properly resamples ALL time-varying data, not just snapshots
    if (
        msv_resolution
        and isinstance(msv_resolution, str)
        and "h" in msv_resolution.lower()
    ):
        n = resample_network_temporal(n, msv_resolution)

    # Prepare network for dispatch (no investment)
    planning_horizons = snakemake.wildcards.get("planning_horizons", None)
    prepare_network(
        n,
        solve_opts=solving.get("options", {}),
        foresight="perfect",  # Always perfect foresight for MSV extraction
        renewable_carriers=[],
        planning_horizons=planning_horizons,
        co2_sequestration_potential=None,
        limit_max_growth=None,
    )

    # Attach config to network for extra_functionality
    n.config = snakemake.config

    # Solve with memory logging
    np.random.seed(solving.get("options", {}).get("seed", 123))

    logging_frequency = solving.get("mem_logging_frequency", 30)
    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=logging_frequency
    ) as mem:
        status, condition = solve_for_msv(
            n,
            solver_name=solver_name,
            solver_options=solver_options,
            io_api=io_api,
            extra_functionality=extra_functionality_msv,
        )

    logger.info(f"Maximum memory usage: {mem.mem_usage}")

    if status != "ok":
        raise RuntimeError(f"MSV extraction failed: {status}, {condition}")

    # Validate MSV and log statistics
    # Get seasonal carriers and threshold from config/params
    seasonal_carriers = snakemake.params.get("seasonal_carriers", [])
    extreme_threshold = snakemake.params.get(
        "msv_extreme_threshold", DEFAULT_MSV_EXTREME_THRESHOLD
    )
    _ = validate_msv(n, seasonal_carriers, extreme_threshold)

    # Verify we have the required dual variables
    if n.stores_t.mu_energy_balance.empty:
        logger.warning(
            "No mu_energy_balance found in stores. "
            "MSV will not be available for rolling horizon optimization."
        )
    else:
        logger.info(
            f"MSV extracted for {len(n.stores_t.mu_energy_balance.columns)} stores"
        )

    # Save network with MSV (dual variables)
    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output.network)
    logger.info(f"Saved network with MSV to {snakemake.output.network}")
