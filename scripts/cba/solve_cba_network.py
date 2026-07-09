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

The optimization is based on the `network.optimize_with_rolling_horizon` method.
Additionally, some extra constraints specified in `solve_network` are added, if
they apply to the dispatch.
"""

import copy
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
from snakemake.utils import update_config
from tqdm.auto import tqdm

from scripts._benchmark import memory_logger
from scripts._helpers import (
    configure_logging,
    get_version,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.solve_network import (
    add_co2_atmosphere_constraint,
    add_import_limit_constraint,
    add_operational_reserve_margin,
    check_objective_value,
    collect_kwargs,
    constrain_dsr_daily_dispatch,
    prepare_network,
)

logger = logging.getLogger(__name__)


def get_components_with_volume_limits(
    n: pypsa.Network,
    type: str,
    carriers: list[str],
) -> pd.Index:
    """
    Return components that have volume limits, given a list of carriers.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network
    type : str
        Component type, e.g. "Generator" or "Link"
    carriers : list of str
        List of carriers to filter components by

    Returns
    -------
    pd.Index
        Index of components that have volume limits and match the given carriers
    """
    static = n.c[type].static
    if "has_volume_limit" not in static.columns:
        return pd.Index([])

    return static.index[static["carrier"].isin(carriers) & static["has_volume_limit"]]


def dispose_model(n: pypsa.Network) -> None:
    """
    Explicitly dispose of Model object.

    This function is relevant when using a single-use license solver for the CBA rolling horizon optimization.
    Without this function, what happens is:
    - First solve (project 1's first rolling horizon): Model object is created and holds the license.
    - After solve completes, Model will be garbage-collected, but license may still be held until next solve.
    - Second solve (project 1's second rolling horizon, or project 2's first rolling horizon): Solver tries
     to create a new Model object, but the license server denies it because the previous object is still active.

    This situation results in the second solve failing with a single-use license error, causing the workflow to fail.

    What this function does is explicitly dispose of the Model after each solve, which releases the license
    and allows subsequent solves to create new Model objects without issue.

    This function should only be called after:
    - A successful solve, OR
    - Computing and printing infeasibilities for a failed solve
    """
    try:
        if (
            n.model is not None
            and hasattr(n.model, "solver_model")
            and n.model.solver_model is not None
        ):
            n.model.solver_model = None
    except Exception as e:
        logger.warning(f"Failed to dispose Model object: {e}")


def extra_functionality(
    n: pypsa.Network,
    snapshots: pd.DatetimeIndex,
    planning_horizons: str | None = None,
) -> None:
    """
    Add custom constraints and functionality for operations network.

    Collects supplementary constraints which will be passed to
    `pypsa.optimization.optimize`.

    If you want to enforce additional custom constraints, this is a good
    location to add them. The arguments `opts` and
    `snakemake.config` are expected to be attached to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance with config and params attributes.
    snapshots : pd.DatetimeIndex
        Simulation timesteps.
    planning_horizons : str or None, optional
        The current planning horizon year or None in perfect foresight.
    """
    config = n.config

    reserve = config["electricity"].get("operational_reserve", {})
    if reserve.get("activate"):
        add_operational_reserve_margin(n, snapshots, config)

    add_co2_atmosphere_constraint(n, snapshots)

    if config["sector"]["imports"]["enable"]:
        add_import_limit_constraint(n, snapshots)

    if config["cba"].get("constrain_dsr", False):
        constrain_dsr_daily_dispatch(n, snapshots)

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
        The PyPSA network instance to optimize.
    snapshots : Sequence or None, optional
        Set of snapshots to consider in the optimization. Default is None.
    horizon : int
        Number of snapshots to consider in each iteration. Default is 100.
    overlap : int
        Number of snapshots to overlap between two iterations. Default is 0.
    **kwargs
        Keyword arguments used by `linopy.Model.solve`, such as `solver_name`.

    Returns
    -------
    tuple[str, str]
        Tuple of (status, condition) from the final optimization window.
    """
    if snapshots is None:
        snapshots: Sequence = n.snapshots

    if horizon <= overlap:
        raise ValueError("overlap must be smaller than horizon")

    assert len(snapshots), "Need at least one snapshot to optimize"

    fallback_solver = kwargs.pop("fallback_solver", None)
    biomass_biogas_slack = n.config["cba"].get("biomass_biogas_slack", 0.4)

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

        # Set per-window energy budgets for biomass/biogas components
        # based on PF dispatch stored in generators_t.p
        for c_name in ["Generator", "Link"]:
            c = n.c[c_name]
            vol_idx = get_components_with_volume_limits(
                n, c_name, ["solid biomass", "biogas"]
            )
            if vol_idx.empty:
                continue
            p_col = "p" if c_name == "Generator" else "p0"
            pf_p = c.dynamic[p_col]
            for comp in vol_idx:
                if comp not in pf_p.columns:
                    continue
                window_energy = pf_p.loc[sns, comp].sum()
                c.static.loc[comp, "e_sum_min"] = (
                    1 - biomass_biogas_slack
                ) * window_energy
                c.static.loc[comp, "e_sum_max"] = window_energy

        status, condition = n.optimize(sns, **kwargs)  # type: ignore

        # If solve is successful, dispose of Model object to release license before next rolling horizon
        if status == "ok":
            dispose_model(n)

        # If solve failed, hold on to the Model object until after IIS is computed in solve_network()
        if status != "ok":
            logger.warning(
                f"Optimization failed with status {status} and condition {condition}"
            )
            # Retry with fallback solver if configured
            if fallback_solver:
                # If solve failed and fallback is configured, dispose of Model before creating a new one.
                dispose_model(n)

                logger.info(
                    f"Retrying window {i + 1}/{len(starting_points)} "
                    f"with fallback solver '{fallback_solver['name']}'"
                )
                retry_kwargs = {**kwargs}
                retry_kwargs["solver_name"] = fallback_solver["name"]
                retry_kwargs["solver_options"] = fallback_solver.get("options", {})
                status, condition = n.optimize(sns, **retry_kwargs)  # type: ignore

                # Only dispose after fallback if it succeeded
                if status == "ok":
                    dispose_model(n)

            if status != "ok":
                logger.warning(f"Fallback also failed: {status} / {condition}")
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
        The PyPSA network instance.
    config : dict
        Configuration dictionary containing solver settings.
    params : dict
        Dictionary of solving parameters.
    solving : dict
        Dictionary of solving options and configuration.
    planning_horizons : str or None, optional
        The current planning horizon year or None in perfect foresight.
    **kwargs
        Additional keyword arguments passed to the solver.

    Raises
    ------
    RuntimeError
        If solving status is infeasible or warning
    ObjectiveValueError
        If objective value differs from expected value
    """
    all_kwargs, _ = collect_kwargs(
        config,
        solving,
        planning_horizons,
        log_fn=kwargs.get("log_fn"),
        mode="rolling_horizon",
    )
    all_kwargs["extra_functionality"] = partial(
        extra_functionality,
        planning_horizons=planning_horizons,
    )

    # Values for horizon and overlap are set in solve_network.collect_kwargs (mode == "rolling_horizon")
    # Thus, we need to override them here with values from the config
    all_kwargs["horizon"] = solving.get("horizon", 168)
    all_kwargs["overlap"] = solving.get("overlap", 1)

    # Configure fallback solver
    fallback_solver = solving.get("fallback_solver", None)
    if fallback_solver:
        fb_options_key = fallback_solver.get("options", "")
        all_kwargs["fallback_solver"] = {
            "name": fallback_solver["name"],
            "options": solving.get("solver_options", {}).get(fb_options_key, {}),
        }

    if all_kwargs.get("solver_name") == "gurobi":
        logging.getLogger("gurobipy").setLevel(logging.CRITICAL)

    # add to network for extra_functionality
    n.config = config
    n.params = params

    status, condition = optimize_with_rolling_horizon(n, **all_kwargs)

    if status != "ok":
        logger.warning(
            f"Solving status '{status}' with termination condition '{condition}'"
        )
        # If infeasible and using Gurobi or Xpress, compute and log infeasibilities before raising error
        if "infeasible" in condition:
            solver_name = solving["solver"]["name"]
            if solver_name in ["gurobi", "xpress"]:
                labels = n.model.compute_infeasibilities()
                logger.info(f"Labels:\n{labels}")
                n.model.print_infeasibilities()
                raise RuntimeError(
                    "Solving status 'infeasible'. Infeasibilities computed."
                )
        raise RuntimeError(
            f"Solving status '{status}' with termination condition '{condition}'."
        )

    check_objective_value(n, solving)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_cba_network",
            run="NT",
            cba_project="t16",
            planning_horizons="2030",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    solving = copy.deepcopy(snakemake.params.solving)
    update_config(solving, snakemake.params.cba_solving)

    np.random.seed(solving["options"].get("seed", 123))

    n = pypsa.Network(snakemake.input.network)
    planning_horizons = snakemake.wildcards.get("planning_horizons", None)

    prepare_network(
        n,
        solve_opts=solving["options"],
        foresight=snakemake.params.foresight,
        renewable_carriers=[],
        planning_horizons=planning_horizons,
        co2_sequestration_potential=None,
        limit_max_growth=None,
        config=snakemake.config,
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

    # Assign meta data to network
    n.meta = dict(
        snakemake.config,
        **dict(wildcards=dict(snakemake.wildcards)),
        version_commit=get_version(),
    )
    n.export_to_netcdf(snakemake.output.network)
