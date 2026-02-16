# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Extract Marginal Storage Values (MSV) via perfect foresight optimization.

Solves the reference network with perfect foresight (full year) to extract
shadow prices from store energy balance constraints. These MSV values capture
the future value of stored energy and guide seasonal storage dispatch in
rolling horizon optimization.

**Inputs**

- ``resources/cba/networks/reference_{planning_horizons}.nc``: Reference network
- ``resources/cba/msv_snapshot_weightings_{planning_horizons}.csv``: Snapshot weightings (optional)

**Outputs**

- ``resources/cba/networks/msv_{planning_horizons}.nc``: Network with MSV in stores_t.mu_energy_balance
"""

import logging

import pypsa
from snakemake.utils import update_config

from scripts._helpers import configure_logging, set_scenario_config
from scripts.prepare_sector_network import set_temporal_aggregation
from scripts.solve_network import prepare_network

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

    # Load base network
    n = pypsa.Network(snakemake.input.network)

    # Optional: resample to coarser resolution for faster solve
    msv_resolution = snakemake.params.get("msv_resolution", False)
    snapshot_weightings = snakemake.input.get("snapshot_weightings", None)
    if msv_resolution:
        n = set_temporal_aggregation(n, msv_resolution, snapshot_weightings)

    # Merge base solving with MSV-specific overrides
    solving = snakemake.params.get("solving", {})
    msv_solving = snakemake.params.get("msv_solving", {})
    update_config(solving, msv_solving)

    solver_name = solving.get("solver", {}).get("name", "highs")
    solver_options_key = solving.get("solver", {}).get("options", "highs-default")

    # Prepare network (e.g., load_shedding setup)
    solve_opts = solving.get("options", {})
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
        logger.error(f"MSV extraction solve failed: {termination_condition}")
        raise RuntimeError(f"MSV extraction solve failed: {termination_condition}")

    logger.info(f"MSV extraction solve completed: {termination_condition}")

    # Log MSV statistics for configured carriers (Store and StorageUnit)
    seasonal_carriers = snakemake.params.get("seasonal_carriers", [])
    accumulator_carriers = snakemake.params.get("accumulator_carriers", [])
    all_msv_carriers = seasonal_carriers + accumulator_carriers

    for component, dual_attr in [
        ("stores", "mu_energy_balance"),
        ("storage_units", "mu_energy_balance"),
    ]:
        dual_df = getattr(getattr(n, f"{component}_t"), dual_attr, None)
        if dual_df is None or dual_df.empty:
            continue
        static = getattr(n, component)
        for carrier in all_msv_carriers:
            idx = static[static.carrier == carrier].index
            idx = idx[idx.isin(dual_df.columns)]
            if idx.empty:
                continue
            msv = dual_df[idx]
            label = "Store" if component == "stores" else "StorageUnit"
            logger.info(
                f"MSV for '{carrier}' ({label}): "
                f"min={msv.min().min():.2f}, max={msv.max().max():.2f}, "
                f"mean={msv.mean().mean():.2f} EUR/MWh"
            )

    # Save network with MSV
    n.export_to_netcdf(snakemake.output.network)
    logger.info(f"Saved MSV network to {snakemake.output.network}")
