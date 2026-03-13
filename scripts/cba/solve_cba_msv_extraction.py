# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Extract marginal storage values via perfect foresight optimization.

Solves the reference network with perfect foresight (full year) to extract
shadow prices from store energy balance constraints. These values capture
the future value of stored energy and guide seasonal storage dispatch in
rolling horizon optimization.

**Inputs**

- ``resources/cba/networks/reference_{planning_horizons}.nc``: Reference network
- ``resources/cba/msv_snapshot_weightings_{planning_horizons}.csv``: Snapshot weightings (optional)

**Outputs**

- ``resources/cba/networks/msv_{planning_horizons}.nc``: Network with marginal storage values in stores_t.mu_energy_balance
"""

import logging

import pypsa
from snakemake.utils import update_config

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)
from scripts.solve_network import prepare_network
from scripts.temporal_aggregation import set_temporal_aggregation

logger = logging.getLogger(__name__)


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

    # Load base network
    n = pypsa.Network(snakemake.input.network)

    # Optional: resample to coarser resolution for faster solve
    msv_resolution = snakemake.params.get("msv_resolution", False)
    snapshot_weightings = snakemake.input.get("snapshot_weightings", None)
    if msv_resolution:
        n = set_temporal_aggregation(n, msv_resolution, snapshot_weightings)

    # Merge CBA-specific solving overrides into the global solving config
    solving = snakemake.params.get("solving", {})
    update_config(solving, snakemake.params.get("cba_solving", {}))

    solver_name = solving.get("solver", {}).get("name", "highs")
    solver_options_key = solving.get("solver", {}).get("options", "highs-default")

    # Prepare network (e.g., load shedding setup)
    solve_opts = solving.get("options", {})
    # Normalize load_shedding: allow bool or dict
    if isinstance(solve_opts.get("load_shedding"), bool):
        solve_opts["load_shedding"] = {"enable": solve_opts["load_shedding"]}
    prepare_network(
        n,
        solve_opts=solve_opts,
        foresight="perfect",
        renewable_carriers=[],
        planning_horizons=snakemake.wildcards.get("planning_horizons", None),
        co2_sequestration_potential=None,
        limit_max_growth=None,
    )

    # Solve with perfect foresight (full year, single optimization)
    # assign_all_duals=True ensures we get mu_energy_balance
    status, termination_condition = n.optimize(
        solver_name=solver_name,
        solver_options=solving.get("solver_options", {}).get(solver_options_key, {}),
        assign_all_duals=True,
    )

    if status != "ok":
        logger.error(f"Extraction solve failed: {termination_condition}")
        # # if the solver is gurobi, print infeasibilities using n.model.print_infeasibilities()
        # if solving.get("solver", {}).get("name", "") == "gurobi":
        #     n.model.print_infeasibilities()
        raise RuntimeError(f"Extraction solve failed: {termination_condition}")

    logger.info(f"Extraction solve completed: {termination_condition}")

    # Save network with marginal storage values
    n.export_to_netcdf(snakemake.output.network)
