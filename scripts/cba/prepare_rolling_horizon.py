# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Prepare network for CBA rolling horizon dispatch.

Takes the reference network and applies modifications for rolling horizon:
- Disable cyclicity for seasonal stores (they receive MSV instead)
- Keep cyclicity for short-term storage (batteries)
- Apply MSV as marginal costs to seasonal stores
- Remove global constraints not needed for CBA
"""

import logging

import pandas as pd
import pypsa
from numpy import inf, isfinite

from scripts._helpers import configure_logging, set_scenario_config
from scripts.cba._helpers import summarize_counts

logger = logging.getLogger(__name__)


def disable_global_constraints(n: pypsa.Network):
    """
    Remove global constraints not needed for CBA rolling horizon.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify
    """
    if "co2_sequestration_limit" in n.global_constraints.index:
        logger.info('Removing GlobalConstraint "co2_sequestration_limit"')
        n.remove("GlobalConstraint", "co2_sequestration_limit")

    if "unsustainable biomass limit" in n.global_constraints.index:
        logger.info('Removing GlobalConstraint "unsustainable biomass limit"')
        n.remove("GlobalConstraint", "unsustainable biomass limit")


def disable_store_cyclicity(
    n: pypsa.Network,
    cyclic_carriers: list[str] | None = None,
    seasonal_carriers: list[str] | None = None,
):
    """
    Disable cyclic constraints for stores and storage units.

    For rolling horizon with MSV, seasonal storage should be non-cyclic
    while short-term storage (batteries) remains cyclic.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify
    cyclic_carriers : list[str], optional
        Carriers that should remain cyclic (e.g., ["battery", "home battery"]).
    seasonal_carriers : list[str], optional
        Carriers that are seasonal (e.g., ["H2 Store", "gas"]). Used for logging.
    """
    if cyclic_carriers is None:
        cyclic_carriers = []
    if seasonal_carriers is None:
        seasonal_carriers = []

    # Disable cyclicity for stores (except cyclic_carriers)
    has_e_cyclic = n.stores["e_cyclic"]
    is_cyclic_carrier = n.stores["carrier"].isin(cyclic_carriers)
    to_disable = has_e_cyclic & ~is_cyclic_carrier

    if to_disable.any():
        stats = summarize_counts(n.stores.loc[to_disable, "carrier"])
        logger.info(f"Disabling e_cyclic for stores:\n{stats}")
        n.stores.loc[to_disable, "e_cyclic"] = False

    if is_cyclic_carrier.any() and has_e_cyclic.any():
        kept_cyclic = has_e_cyclic & is_cyclic_carrier
        if kept_cyclic.any():
            stats = summarize_counts(n.stores.loc[kept_cyclic, "carrier"])
            logger.info(f"Keeping e_cyclic=True for short-term storage:\n{stats}")

    # Log seasonal carriers that will receive MSV
    if seasonal_carriers:
        is_seasonal = n.stores["carrier"].isin(seasonal_carriers)
        if is_seasonal.any():
            stats = summarize_counts(n.stores.loc[is_seasonal, "carrier"])
            logger.info(f"Seasonal stores (will receive MSV):\n{stats}")

    has_e_cyclic_per_period = n.stores["e_cyclic_per_period"]
    to_disable_per_period = has_e_cyclic_per_period & ~is_cyclic_carrier

    if to_disable_per_period.any():
        stats = summarize_counts(n.stores.loc[to_disable_per_period, "carrier"])
        logger.info(f"Disabling e_cyclic_per_period for stores:\n{stats}")
        n.stores.loc[to_disable_per_period, "e_cyclic_per_period"] = False

    # Disable cyclicity for storage units (except cyclic_carriers)
    has_cyclic_soc = n.storage_units["cyclic_state_of_charge"]
    is_cyclic_carrier_su = n.storage_units["carrier"].isin(cyclic_carriers)
    to_disable_su = has_cyclic_soc & ~is_cyclic_carrier_su

    if to_disable_su.any():
        stats = summarize_counts(n.storage_units.loc[to_disable_su, "carrier"])
        logger.info(f"Disabling cyclic_state_of_charge for storage units:\n{stats}")
        n.storage_units.loc[to_disable_su, "cyclic_state_of_charge"] = False

    if is_cyclic_carrier_su.any() and has_cyclic_soc.any():
        kept_cyclic_su = has_cyclic_soc & is_cyclic_carrier_su
        if kept_cyclic_su.any():
            stats = summarize_counts(n.storage_units.loc[kept_cyclic_su, "carrier"])
            logger.info(
                f"Keeping cyclic_state_of_charge=True for short-term storage:\n{stats}"
            )


def resample_msv_to_target(
    msv: pd.DataFrame,
    target_snapshots: pd.DatetimeIndex,
    method: str = "ffill",
) -> pd.DataFrame:
    """
    Resample MSV from extraction resolution to target network resolution.

    Parameters
    ----------
    msv : pd.DataFrame
        MSV data from extraction (e.g., 24H resolution)
    target_snapshots : pd.DatetimeIndex
        Target snapshots (e.g., 3H resolution)
    method : str, optional
        Resampling method:
        - "ffill": Forward fill - each MSV value applies until the next one
        - "interpolate": Linear interpolation between MSV values

    Returns
    -------
    pd.DataFrame
        MSV resampled to target resolution
    """
    if method == "interpolate":
        # Combine indices and interpolate
        combined_index = msv.index.union(target_snapshots).sort_values()
        msv_resampled = msv.reindex(combined_index).interpolate(method="time")
        msv_resampled = msv_resampled.reindex(target_snapshots)
    else:
        # Default: forward fill
        if method != "ffill":
            logger.warning(f"Unknown resample method '{method}', using 'ffill'")
        msv_resampled = msv.reindex(target_snapshots, method="ffill")

    # Fill any remaining NaNs at the start with backward fill
    msv_resampled = msv_resampled.bfill()

    return msv_resampled


def disable_volume_limits(n: pypsa.Network):
    """
    Disable minimum energy production limits (e_sum_min) for generators and links.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify
    """
    for c in n.components[{"Generator", "Link"}]:
        has_e_sum_min = isfinite(c.static.get("e_sum_min", []))
        if has_e_sum_min.any():
            stats = summarize_counts(c.static.loc[has_e_sum_min, "carrier"])
            logger.info(f"Disabling e_sum_min volume limits of:\n{stats}")
            c.static.loc[has_e_sum_min, "e_sum_min"] = -inf


def apply_msv_to_network(
    n: pypsa.Network,
    n_msv: pypsa.Network,
    seasonal_carriers: list[str],
    resample_method: str = "ffill",
) -> None:
    """
    Apply Marginal Storage Values (MSV) to seasonal stores.

    The MSV (mu_energy_balance dual variable) from perfect foresight is applied
    as marginal_cost to stores, guiding rolling horizon dispatch to account for
    seasonal storage value.

    Parameters
    ----------
    n : pypsa.Network
        Target network for rolling horizon optimization
    n_msv : pypsa.Network
        Network with MSV from perfect foresight solve
    seasonal_carriers : list[str]
        List of carrier names that should receive MSV
    resample_method : str, optional
        Method for resampling MSV to target resolution ("ffill" or "interpolate")
    """
    if not seasonal_carriers:
        logger.info("No seasonal carriers specified, skipping MSV application")
        return

    if n_msv.stores_t.mu_energy_balance.empty:
        logger.warning(
            "MSV network has no mu_energy_balance. "
            "Ensure MSV extraction was solved with assign_all_duals=True. "
            "Skipping MSV application."
        )
        return

    # Check if resampling is needed
    needs_resample = not n.snapshots.equals(n_msv.snapshots)
    if needs_resample:
        logger.info(
            f"MSV resampling: {len(n_msv.snapshots)} → {len(n.snapshots)} snapshots "
            f"(method: {resample_method})"
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
            logger.warning(f"No MSV data for stores with carrier '{carrier}'")
            continue

        # Extract MSV
        msv = n_msv.stores_t.mu_energy_balance[msv_stores_i].copy()

        # Resample if needed
        if needs_resample:
            msv = resample_msv_to_target(msv, n.snapshots, method=resample_method)

        # Apply MSV as marginal_cost
        if n.stores_t.marginal_cost.empty:
            n.stores_t.marginal_cost = msv
        else:
            for store in msv.columns:
                n.stores_t.marginal_cost[store] = msv[store]

        logger.info(
            f"Applied MSV to {len(msv_stores_i)} stores with carrier '{carrier}'. "
            f"MSV range: [{msv.min().min():.2f}, {msv.max().max():.2f}] EUR/MWh"
        )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_rolling_horizon",
            planning_horizons="2030",
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the CBA reference network (already has fixed capacities and hurdle costs)
    n = pypsa.Network(snakemake.input.get("network", None))

    # Get storage carrier settings from config
    cyclic_carriers = snakemake.params.get("cyclic_carriers", [])
    seasonal_carriers = snakemake.params.get("seasonal_carriers", [])
    disable_store_cyclicity(
        n, cyclic_carriers=cyclic_carriers, seasonal_carriers=seasonal_carriers
    )

    disable_global_constraints(n)

    # disable volume limits
    disable_volume_limits(n)

    # Load MSV network and apply MSV to seasonal stores
    msv_network_path = snakemake.input.get("network_msv", None)
    if msv_network_path:
        n_msv = pypsa.Network(msv_network_path)
        resample_method = snakemake.params.get("msv_resample_method", "ffill")
        apply_msv_to_network(n, n_msv, seasonal_carriers, resample_method)
    else:
        logger.info("No MSV network provided, skipping MSV application")

    # Save simplified network
    n.export_to_netcdf(snakemake.output.network)
