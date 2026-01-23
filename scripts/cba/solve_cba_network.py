# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Solves optimal operation in rolling horizons for fixed capacities.

This script is used for optimizing the electrical network as well as the
sector coupled network.

Description
-----------

The optimization is based on the :func:`network.optimize_with_rolling_horizon` method.
Additionally, some extra constraints specified in :mod:`solve_network` are added, if
they apply to the dispatch.
"""

import importlib
import logging
import os
import sys
from collections.abc import Sequence
from functools import partial
from typing import Any

import numpy as np
import pandas as pd
import pypsa
from linopy.remote.oetc import OetcCredentials, OetcHandler, OetcSettings
from snakemake.utils import update_config
from tqdm.auto import tqdm

from scripts._benchmark import memory_logger
from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.solve_network import (
    add_co2_atmosphere_constraint,
    add_import_limit_constraint,
    add_operational_reserve_margin,
    check_objective_value,
    prepare_network,
)

logger = logging.getLogger(__name__)

# Default threshold for warning about extreme MSV values (EUR/MWh)
DEFAULT_MSV_EXTREME_THRESHOLD = 1000.0


def resample_msv_to_target(
    msv: pd.DataFrame,
    target_snapshots: pd.DatetimeIndex,
    method: str = "ffill",
    rolling_window: int | None = None,
) -> pd.DataFrame:
    """
    Resample MSV from one resolution to target resolution.

    Parameters
    ----------
    msv : pd.DataFrame
        MSV data (e.g., 3-hourly from MSV extraction)
    target_snapshots : pd.DatetimeIndex
        Target snapshots (e.g., hourly for rolling horizon)
    method : str, optional
        Resampling method. Options:
        - "ffill" (default): Forward fill - each MSV value applies until the next one.
          Appropriate for shadow prices which represent discrete time periods.
        - "bfill": Backward fill - each MSV value applies from the previous one.
        - "interpolate": Linear interpolation - smoother but may not reflect discrete
          price steps.
        - "nearest": Use nearest MSV value by time.
        - "rolling_mean": Apply rolling mean smoothing after forward fill.
          Reduces volatility while preserving trends. Requires rolling_window.
    rolling_window : int, optional
        Window size for rolling mean (number of snapshots). Only used when
        method="rolling_mean". Default is None (must be specified if using rolling_mean).

    Returns
    -------
    pd.DataFrame
        MSV resampled to target resolution

    Notes
    -----
    Handles both upsampling (MSV coarser than target, typical case) and
    downsampling (MSV finer than target, unusual but possible).
    """
    logger.info(f"Resampling MSV using method='{method}'")

    # Check if target is coarser than MSV (downsampling case)
    if len(target_snapshots) < len(msv):
        logger.info(
            f"Downsampling MSV from {len(msv)} to {len(target_snapshots)} snapshots "
            f"(using mean aggregation)"
        )
        # For downsampling, use mean to aggregate MSV values
        msv_resampled = msv.reindex(target_snapshots, method="nearest")
        return msv_resampled

    # Upsampling case (typical: MSV at 3H, target at 1H)
    if method == "interpolate":
        # Linear interpolation - smoother but may not reflect discrete price steps
        combined_index = msv.index.union(target_snapshots).sort_values()
        msv_resampled = msv.reindex(combined_index).interpolate(method="time")
        msv_resampled = msv_resampled.reindex(target_snapshots)

    elif method == "bfill":
        # Backward fill - each MSV value applies from the previous one
        msv_resampled = msv.reindex(target_snapshots, method="bfill")

    elif method == "nearest":
        # Use nearest MSV value by time
        msv_resampled = msv.reindex(target_snapshots, method="nearest")

    elif method == "rolling_mean":
        # Forward fill first, then apply rolling mean to smooth
        if rolling_window is None:
            raise ValueError("rolling_window must be specified when method='rolling_mean'")
        msv_resampled = msv.reindex(target_snapshots, method="ffill")
        msv_resampled = msv_resampled.rolling(window=rolling_window, min_periods=1, center=True).mean()
        logger.info(f"Applied rolling mean with window={rolling_window} snapshots")

    else:
        # Default: Forward fill - each MSV value applies until the next one
        # This is most appropriate for shadow prices which represent discrete time periods
        if method != "ffill":
            logger.warning(f"Unknown resampling method '{method}', falling back to 'ffill'")
        msv_resampled = msv.reindex(target_snapshots, method="ffill")

    # Fill any remaining NaNs at the beginning with backward fill
    msv_resampled = msv_resampled.bfill()

    return msv_resampled


def apply_marginal_storage_values(
    n: pypsa.Network,
    n_msv: pypsa.Network,
    seasonal_carriers: list[str],
    msv_extreme_threshold: float = DEFAULT_MSV_EXTREME_THRESHOLD,
    resample_method: str = "ffill",
    rolling_window: int | None = None,
) -> None:
    """
    Apply marginal storage values (MSV) to seasonal stores.

    MSV are dual variables from the store energy balance constraint obtained from
    a full-year perfect foresight optimization. They capture the future value of
    stored energy. By setting these as marginal costs on stores, the rolling horizon
    optimization can make decisions that account for seasonal patterns.

    The MSV is extracted from n_msv.stores_t.mu_energy_balance, which is the dual
    variable of the Store-energy_balance constraint.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify (for rolling horizon optimization)
    n_msv : pypsa.Network
        Network solved with perfect foresight containing mu_energy_balance
    seasonal_carriers : list[str]
        List of store carrier names that should receive MSV
        (e.g., ["H2 Store", "gas", "co2 sequestered"])
    msv_extreme_threshold : float, optional
        Threshold for warning about extreme MSV values (EUR/MWh). Default is 1000.
    resample_method : str, optional
        Method for resampling MSV to target resolution. Options: "ffill" (default),
        "bfill", "interpolate", "nearest", "rolling_mean". Default is "ffill".
    rolling_window : int, optional
        Window size for rolling mean (snapshots). Required if resample_method="rolling_mean".

    Notes
    -----
    If the MSV network has different temporal resolution than the target network,
    the MSV values are resampled using the specified method.

    This implements the concept from Lagrangian relaxation where the MSV acts
    as a penalty/reward for storing energy, internalizing future scarcity even
    when optimizing over limited horizons.

    References
    ----------
    - Tom Brown's optimization lecture: https://nworbmot.org/courses/es-24/es-7-optimisation.pdf
    - PyPSA documentation on storage optimization
    """
    if not seasonal_carriers:
        logger.info("No seasonal carriers specified, skipping MSV assignment")
        return

    # Check if MSV network has mu_energy_balance (the correct dual variable)
    if n_msv.stores_t.mu_energy_balance.empty:
        logger.warning(
            "MSV network has no mu_energy_balance. "
            "Ensure it was solved with assign_all_duals=True. "
            "Skipping MSV assignment."
        )
        return

    logger.info(
        f"Applying MSV from network with {len(n_msv.snapshots)} snapshots "
        f"to network with {len(n.snapshots)} snapshots"
    )

    for carrier in seasonal_carriers:
        # Get stores with this carrier in target network
        stores_i = n.stores[n.stores.carrier == carrier].index

        if stores_i.empty:
            logger.debug(f"No stores found with carrier '{carrier}'")
            continue

        # Find corresponding stores in MSV network
        msv_stores_i = stores_i[stores_i.isin(n_msv.stores_t.mu_energy_balance.columns)]

        if msv_stores_i.empty:
            logger.warning(
                f"No MSV data found for stores with carrier '{carrier}'. "
                "Skipping these stores."
            )
            continue

        # Extract MSV from mu_energy_balance
        msv = n_msv.stores_t.mu_energy_balance[msv_stores_i].copy()

        # Resample if temporal resolutions differ
        if not n.snapshots.equals(n_msv.snapshots):
            logger.info(
                f"Resampling MSV for '{carrier}' from {len(n_msv.snapshots)} "
                f"to {len(n.snapshots)} snapshots"
            )
            msv = resample_msv_to_target(
                msv, n.snapshots, method=resample_method, rolling_window=rolling_window
            )

        # Ensure all target stores are included (some may not be in MSV network)
        # For stores not in MSV network (e.g., new PINT storage projects),
        # we keep their existing marginal_cost (typically 0)
        missing_stores = stores_i.difference(msv_stores_i)
        if not missing_stores.empty:
            logger.warning(
                f"{len(missing_stores)} stores with carrier '{carrier}' not found "
                f"in MSV network (possibly new PINT storage). "
                f"These stores will use marginal_cost=0, which may lead to "
                f"suboptimal seasonal dispatch. Consider running MSV extraction "
                f"with these stores included for better results."
            )

        # Set the marginal cost of stores to the MSV
        # The MSV (mu_energy_balance) represents the marginal value of stored energy.
        # Positive MSV means stored energy has value (future scarcity expected).
        # By setting this as marginal_cost, the optimizer will:
        # - Avoid discharging when MSV is high (energy is valuable)
        # - Prefer charging when MSV is high (building reserves)
        # - Discharge when MSV is low (energy less valuable in future)

        # Validate MSV values before applying
        if (msv.abs() > msv_extreme_threshold).any().any():
            extreme_count = (msv.abs() > msv_extreme_threshold).sum().sum()
            logger.warning(
                f"Extreme MSV values (|MSV| > {msv_extreme_threshold} EUR/MWh) "
                f"detected for '{carrier}': {extreme_count} occurrences. "
                "This may cause unusual dispatch behavior."
            )

        if n.stores_t.marginal_cost.empty:
            n.stores_t.marginal_cost = msv
        else:
            # Update existing marginal_cost DataFrame with MSV values
            for store in msv.columns:
                n.stores_t.marginal_cost[store] = msv[store]

        logger.info(
            f"Applied MSV to {len(msv_stores_i)} stores with carrier '{carrier}'. "
            f"MSV range: [{msv.min().min():.2f}, {msv.max().max():.2f}] EUR/MWh"
        )

    # Set initial energy levels from MSV network for seasonal stores
    # This ensures consistency between the full-year solution and rolling horizon
    for carrier in seasonal_carriers:
        stores_i = n.stores[n.stores.carrier == carrier].index
        if stores_i.empty:
            continue

        # Get initial energy from MSV network's first snapshot
        msv_stores_i = stores_i[stores_i.isin(n_msv.stores.index)]
        if not msv_stores_i.empty and not n_msv.stores_t.e.empty:
            first_snapshot = n_msv.snapshots[0]
            if first_snapshot in n_msv.stores_t.e.index:
                initial_e = n_msv.stores_t.e.loc[first_snapshot, msv_stores_i]
                n.stores.loc[msv_stores_i, "e_initial"] = initial_e
                logger.info(
                    f"Set initial energy for {len(msv_stores_i)} '{carrier}' stores "
                    f"from MSV network (first snapshot: {first_snapshot})"
                )


def extra_functionality(
    n: pypsa.Network,
    snapshots: pd.DatetimeIndex,
    planning_horizons: str | None = None,
) -> None:
    """
    Add custom constraints and functionality for operations network

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance with config and params attributes
    snapshots : pd.DatetimeIndex
        Simulation timesteps
    planning_horizons : str, optional
        The current planning horizon year or None in perfect foresight

    Collects supplementary constraints which will be passed to
    ``pypsa.optimization.optimize``.

    If you want to enforce additional custom constraints, this is a good
    location to add them. The arguments ``opts`` and
    ``snakemake.config`` are expected to be attached to the network.
    """
    config = n.config

    reserve = config["electricity"].get("operational_reserve", {})
    if reserve.get("activate"):
        add_operational_reserve_margin(n, snapshots, config)

    add_co2_atmosphere_constraint(n, snapshots)

    if config["sector"]["imports"]["enable"]:
        add_import_limit_constraint(n, snapshots)

    if n.params.custom_extra_functionality:
        source_path = n.params.custom_extra_functionality
        assert os.path.exists(source_path), f"{source_path} does not exist"
        sys.path.append(os.path.dirname(source_path))
        module_name = os.path.splitext(os.path.basename(source_path))[0]
        module = importlib.import_module(module_name)
        custom_extra_functionality = getattr(module, module_name)
        custom_extra_functionality(n, snapshots, snakemake)  # pylint: disable=E0601


# TODO should be upstreamed back and replace pypsa.optimization.abstract.optimize_with_rolling_horizon, which
# has currently broken status updates.
def optimize_with_rolling_horizon(
    n: pypsa.Network,
    snapshots: Sequence | None = None,
    horizon: int = 100,
    overlap: int = 0,
    **kwargs: Any,
) -> tuple[str, str]:
    """
    Optimizes the network in a rolling horizon fashion.

    Parameters
    ----------
    n : pypsa.Network
    snapshots : list-like
        Set of snapshots to consider in the optimization. The default is None.
    horizon : int
        Number of snapshots to consider in each iteration. Defaults to 100.
    overlap : int
        Number of snapshots to overlap between two iterations. Defaults to 0.
    **kwargs:
        Keyword argument used by `linopy.Model.solve`, such as `solver_name`,

    Returns
    -------
    tuple[str, str]
    """
    if snapshots is None:
        snapshots = n.snapshots

    if horizon <= overlap:
        raise ValueError("overlap must be smaller than horizon")

    assert len(snapshots), "Need at least one snapshot to optimize"

    starting_points = range(0, len(snapshots), horizon - overlap)
    for i, start in tqdm(enumerate(starting_points), total=len(starting_points)):
        end = min(len(snapshots), start + horizon)
        sns = snapshots[start:end]

        msg = f"Optimizing network for snapshot horizon [{sns[0]}:{sns[-1]}] ({i + 1}/{len(starting_points)})."
        logger.info(msg)
        if log_fn := kwargs.get("log_fn"):
            with open(log_fn, "a") as f:
                print(20 * "=", file=f)
                print(msg, file=f)
                print(20 * "=" + "\n", file=f)

        if i:
            if not n.stores.empty:
                n.stores.e_initial = n.stores_t.e.loc[snapshots[start - 1]]
            if not n.storage_units.empty:
                n.storage_units.state_of_charge_initial = (
                    n.storage_units_t.state_of_charge.loc[snapshots[start - 1]]
                )

        status, condition = n.optimize(sns, **kwargs)  # type: ignore
        if status != "ok":
            logger.warning(
                f"Optimization failed with status {status} and condition {condition}"
            )
            return status, condition

    return status, condition  # pyright: ignore[reportPossiblyUnboundVariable]


def solve_network(
    n: pypsa.Network,
    config: dict,
    params: dict,
    solving: dict,
    planning_horizons: str | None = None,
    **kwargs,
) -> None:
    """
    Solve network optimization problem.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    config : Dict
        Configuration dictionary containing solver settings
    params : Dict
        Dictionary of solving parameters
    solving : Dict
        Dictionary of solving options and configuration
    rule_name : str, optional
        Name of the snakemake rule being executed
    planning_horizons : str, optional
        The current planning horizon year or None in perfect foresight
    **kwargs
        Additional keyword arguments passed to the solver

    Returns
    -------
    n : pypsa.Network
        Solved network instance
    status : str
        Solution status
    condition : str
        Termination condition

    Raises
    ------
    RuntimeError
        If solving status is infeasible or warning
    ObjectiveValueError
        If objective value differs from expected value
    """
    set_of_options = solving["solver"]["options"]
    cf_solving = solving["options"]

    kwargs["multi_investment_periods"] = config["foresight"] == "perfect"
    kwargs["solver_options"] = (
        solving["solver_options"][set_of_options] if set_of_options else {}
    )
    kwargs["solver_name"] = solving["solver"]["name"]
    kwargs["extra_functionality"] = partial(
        extra_functionality,
        planning_horizons=planning_horizons,
    )
    kwargs["transmission_losses"] = cf_solving.get("transmission_losses", False)
    kwargs["linearized_unit_commitment"] = cf_solving.get(
        "linearized_unit_commitment", False
    )
    kwargs["assign_all_duals"] = cf_solving.get("assign_all_duals", False)
    kwargs["io_api"] = cf_solving.get("io_api", None)

    oetc = solving.get("oetc", None)
    if oetc:
        oetc["credentials"] = OetcCredentials(
            email=os.environ["OETC_EMAIL"], password=os.environ["OETC_PASSWORD"]
        )
        oetc["solver"] = kwargs["solver_name"]
        oetc["solver_options"] = kwargs["solver_options"]
        oetc_settings = OetcSettings(**oetc)
        oetc_handler = OetcHandler(oetc_settings)
        kwargs["remote"] = oetc_handler

    kwargs["model_kwargs"] = cf_solving.get("model_kwargs", {})
    kwargs["keep_files"] = cf_solving.get("keep_files", False)

    if kwargs["solver_name"] == "gurobi":
        logging.getLogger("gurobipy").setLevel(logging.CRITICAL)

    # add to network for extra_functionality
    n.config = config
    n.params = params

    kwargs["horizon"] = cf_solving.get("horizon", 24 * 7)
    kwargs["overlap"] = cf_solving.get("overlap", 0)

    status, condition = optimize_with_rolling_horizon(n, **kwargs)

    if status != "ok":
        logger.warning(
            f"Solving status '{status}' with termination condition '{condition}'"
        )
    check_objective_value(n, solving)

    if "warning" in condition:
        raise RuntimeError("Solving status 'warning'. Discarding solution.")

    if "infeasible" in condition:
        labels = n.model.compute_infeasibilities()
        logger.info(f"Labels:\n{labels}")
        n.model.print_infeasibilities()
        raise RuntimeError("Solving status 'infeasible'. Infeasibilities computed.")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_cba_network",
            cba_method="toot",
            name="reference",
            planning_horizons="2030",
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    solving = snakemake.params.solving
    update_config(solving, snakemake.params.cba_solving)

    np.random.seed(solving["options"].get("seed", 123))

    n = pypsa.Network(snakemake.input.network)
    planning_horizons = snakemake.wildcards.get("planning_horizons", None)

    # Load MSV network (solved with perfect foresight to extract marginal storage values)
    n_msv = pypsa.Network(snakemake.input.msv_network)

    # Apply marginal storage values (MSV) to seasonal stores
    seasonal_carriers = snakemake.params.get("seasonal_carriers", [])
    msv_extreme_threshold = snakemake.params.get(
        "msv_extreme_threshold", DEFAULT_MSV_EXTREME_THRESHOLD
    )
    resample_method = snakemake.params.get("msv_resample_method", "ffill")
    rolling_window = snakemake.params.get("msv_rolling_window", None)
    if seasonal_carriers:
        apply_marginal_storage_values(
            n,
            n_msv,
            seasonal_carriers,
            msv_extreme_threshold,
            resample_method=resample_method,
            rolling_window=rolling_window,
        )

    prepare_network(
        n,
        solve_opts=solving["options"],
        foresight=snakemake.params.foresight,
        renewable_carriers=[],
        planning_horizons=planning_horizons,
        co2_sequestration_potential=None,
        limit_max_growth=None,
    )

    logging_frequency = solving.get("mem_logging_frequency", 30)
    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=logging_frequency
    ) as mem:
        solve_network(
            n,
            config=snakemake.config,
            params=snakemake.params,
            solving=solving,
            planning_horizons=planning_horizons,
            log_fn=snakemake.log.solver,
        )

    logger.info(f"Maximum memory usage: {mem.mem_usage}")

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output.network)
