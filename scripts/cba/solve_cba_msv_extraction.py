# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Extract Marginal Storage Values (MSV) via perfect foresight optimization.

This script solves the CBA base network with perfect foresight (full year) to extract
the dual variables (shadow prices) from the store energy balance constraint. These MSV
values capture the future value of stored energy and are used to guide seasonal storage
dispatch in rolling horizon optimization.

The MSV is extracted from n.stores_t.mu_energy_balance, which is the dual variable
of the Store-energy_balance constraint.

**Key Features**

- Optional temporal resampling before solve (e.g., 3H → 24H for faster solve)
- Extracts mu_energy_balance for seasonal stores
- Stores with cyclic constraint enabled to get meaningful shadow prices

**Inputs**

- ``resources/cba/networks/base_{planning_horizons}.nc``: CBA base network

**Outputs**

- ``resources/cba/networks/msv_{planning_horizons}.nc``: Network with MSV in stores_t.mu_energy_balance
"""

import logging

import pypsa

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_cba_msv_extraction",
            planning_horizons="2030",
            run="test-sector-tyndp",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load base network
    n = pypsa.Network(snakemake.input.network)
    logger.info(f"Loaded base network with {len(n.snapshots)} snapshots")

    # Optional: resample to coarser resolution for faster solve
    msv_resolution = snakemake.params.get("msv_resolution", False)
    if msv_resolution:
        n = resample_network(n, msv_resolution)

    # Get solver settings
    solving = snakemake.params.get("solving", {})
    solver_name = solving.get("solver", {}).get("name", "highs")
    solver_options = solving.get("solver", {}).get("options", "highs-default")

    logger.info(f"Solving MSV extraction with {solver_name} ({solver_options})")

    # Solve with perfect foresight (full year, single optimization)
    # assign_all_duals=True ensures we get mu_energy_balance
    status, termination_condition = n.optimize(
        solver_name=solver_name,
        solver_options=solving.get("solver_options", {}).get(solver_options, {}),
        assign_all_duals=True,
    )

    if status != "ok":
        logger.error(f"MSV extraction solve failed: {termination_condition}")
        raise RuntimeError(f"MSV extraction solve failed: {termination_condition}")

    logger.info(f"MSV extraction solve completed: {termination_condition}")

    # Log MSV statistics for seasonal carriers
    seasonal_carriers = snakemake.params.get("seasonal_carriers", [])
    if not n.stores_t.mu_energy_balance.empty:
        for carrier in seasonal_carriers:
            stores_i = n.stores[n.stores.carrier == carrier].index
            stores_i = stores_i[stores_i.isin(n.stores_t.mu_energy_balance.columns)]
            if stores_i.empty:
                continue
            msv = n.stores_t.mu_energy_balance[stores_i]
            logger.info(
                f"MSV for '{carrier}': "
                f"min={msv.min().min():.2f}, max={msv.max().max():.2f}, "
                f"mean={msv.mean().mean():.2f} EUR/MWh"
            )
    else:
        logger.warning("No mu_energy_balance extracted - check solver settings")

    # Save network with MSV
    n.export_to_netcdf(snakemake.output.network)
    logger.info(f"Saved MSV network to {snakemake.output.network}")
