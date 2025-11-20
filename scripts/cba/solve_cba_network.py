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
import yaml
from linopy.remote.oetc import OetcCredentials, OetcHandler, OetcSettings

from scripts._benchmark import memory_logger
from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.solve_network import (
    add_chp_constraints,
    add_co2_atmosphere_constraint,
    add_import_limit_constraint,
    add_operational_reserve_margin,
    check_objective_value,
    prepare_network,
)

logger = logging.getLogger(__name__)


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
    add_chp_constraints(n)

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
        snapshots: Sequence = n.snapshots

    if horizon <= overlap:
        raise ValueError("overlap must be smaller than horizon")

    assert len(snapshots), "Need at least one snapshot to optimize"

    starting_points = range(0, len(snapshots), horizon - overlap)
    for i, start in enumerate(starting_points):
        end = min(len(snapshots), start + horizon)
        sns = snapshots[start:end]
        logger.info(
            f"Optimizing network for snapshot horizon [{sns[0]}:{sns[-1]}] ({i + 1}/{len(starting_points)})."
        )

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
    offshore_zone_trajectories_fn : str, optional
        Path to DataFrame containing the offshore zone potentials trajectories
    renewable_carriers_tyndp : list[str], optional
        List of TYNDP renewable carriers
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

    kwargs["horizon"] = cf_solving.get("horizon", 24 * 7 // params.time_resolution)
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
            configfiles=["config/test/config.tyndp.yaml"],
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    solve_opts = snakemake.params.solving["options"]

    np.random.seed(solve_opts.get("seed", 123))

    n = pypsa.Network(snakemake.input.network)
    planning_horizons = snakemake.wildcards.get("planning_horizons", None)

    # should probably already have happened in simplify_cba_network, applying twice does not hurt, though
    n.optimize.fix_optimal_capacities()

    prepare_network(
        n,
        solve_opts=snakemake.params.solving["options"],
        foresight=snakemake.params.foresight,
        renewable_carriers=snakemake.params.renewable_carriers,
        planning_horizons=planning_horizons,
        co2_sequestration_potential=snakemake.params["co2_sequestration_potential"],
        limit_max_growth=snakemake.params.get("sector", {}).get("limit_max_growth"),
    )

    logging_frequency = snakemake.config.get("solving", {}).get(
        "mem_logging_frequency", 30
    )
    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=logging_frequency
    ) as mem:
        solve_network(
            n,
            config=snakemake.config,
            params=snakemake.params,
            solving=snakemake.params.solving,
            planning_horizons=planning_horizons,
            log_fn=snakemake.log.solver,
        )

    logger.info(f"Maximum memory usage: {mem.mem_usage}")

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output.network)

    with open(snakemake.output.config, "w") as file:
        yaml.dump(
            n.meta,
            file,
            default_flow_style=False,
            allow_unicode=True,
            sort_keys=False,
        )
