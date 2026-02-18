# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Prepare network for CBA rolling horizon dispatch.

Modifications applied:
- cyclic_carriers: keep cyclicity (short-term, e.g., batteries)
- seasonal_carriers: disable cyclicity, apply MSV + initial state from PF
- accumulator_carriers: disable cyclicity, apply MSV, start empty
- Remove global constraints not needed for CBA
"""

import logging

import numpy as np
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
    accumulator_carriers: list[str] | None = None,
):
    """
    Disable cyclic constraints for stores and storage units.

    cyclic_carriers remain cyclic; all others become non-cyclic.
    seasonal_carriers and accumulator_carriers are used for logging only.
    """
    if cyclic_carriers is None:
        cyclic_carriers = []
    if seasonal_carriers is None:
        seasonal_carriers = []
    if accumulator_carriers is None:
        accumulator_carriers = []

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
            # TODO this should be later changed by adding MSV as well
            c.dynamic.p_set.loc[:, c.static[has_e_sum_min].index] = c.dynamic.p.loc[
                :, c.static[has_e_sum_min].index
            ]


def apply_msv_to_network(
    n: pypsa.Network,
    n_msv: pypsa.Network,
    seasonal_carriers: list[str],
    accumulator_carriers: list[str] | None = None,
    resample_method: str = "ffill",
) -> None:
    """
    Apply MSV (mu_energy_balance duals) as marginal_cost to Store and StorageUnit
    components whose carrier is in seasonal_carriers or accumulator_carriers.

    Only explicitly configured carriers receive MSV.
    """
    if accumulator_carriers is None:
        accumulator_carriers = []

    all_msv_carriers = seasonal_carriers + accumulator_carriers

    if not all_msv_carriers:
        logger.info(
            "No seasonal or accumulator carriers specified, skipping MSV application"
        )
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

    for c in ["Store", "StorageUnit"]:
        # Filter to configured carriers only (same logic for both component types)
        s_i = n_msv.c[c].static[n_msv.c[c].static.carrier.isin(all_msv_carriers)].index
        # get shadow prices
        msv = n_msv.c[c].dynamic["mu_energy_balance"][s_i]
        # Resample if needed
        if needs_resample:
            msv = resample_msv_to_target(msv, n.snapshots, method=resample_method)
        # set shadow prices as marginal cost
        n.c[c].dynamic["marginal_cost"].loc[:, s_i] = msv

        logger.info(
            f"Applied MSV to {len(s_i)} {c}. "
            f"MSV range: [{msv.min().min():.2f}, {msv.max().max():.2f}] EUR/MWh"
        )


def set_initial_state_from_pf(
    n: pypsa.Network,
    n_msv: pypsa.Network,
    seasonal_carriers: list[str],
) -> None:
    """
    Set initial storage state from PF solution for seasonal_carriers.

    Store: e_initial = PF e(t=-1).
    StorageUnit: state_of_charge_initial = PF soc(t=-1).

    Not applied to accumulator_carriers (they start empty).
    """
    # Set e_initial for seasonal stores from PF e(t=-1)
    if not n_msv.stores_t.e.empty:
        pf_e_initial = n_msv.stores_t.e.iloc[-1]
        # Filter to seasonal carriers only
        seasonal_stores = n.stores[n.stores.carrier.isin(seasonal_carriers)].index
        common_stores = seasonal_stores.intersection(pf_e_initial.index)
        if len(common_stores) > 0:
            n.stores.loc[common_stores, "e_initial"] = pf_e_initial.loc[common_stores]
            stats = summarize_counts(n.stores.loc[common_stores, "carrier"])
            logger.info(
                f"Set e_initial from PF for {len(common_stores)} stores:\n{stats}"
            )

    # Set state_of_charge_initial for storage units from PF soc(t=-1)
    if not n_msv.storage_units_t.state_of_charge.empty:
        pf_soc_initial = n_msv.storage_units_t.state_of_charge.iloc[-1]
        # Filter to seasonal carriers only
        seasonal_sus = n.storage_units[
            n.storage_units.carrier.isin(seasonal_carriers)
        ].index
        common_sus = seasonal_sus.intersection(pf_soc_initial.index)
        if len(common_sus) > 0:
            n.storage_units.loc[common_sus, "state_of_charge_initial"] = (
                pf_soc_initial.loc[common_sus]
            )
            stats = summarize_counts(n.storage_units.loc[common_sus, "carrier"])
            logger.info(
                f"Set state_of_charge_initial from PF for {len(common_sus)} storage units:\n{stats}"
            )


def fix_reservoir_soc_at_boundaries(
    n: pypsa.Network,
    n_msv: pypsa.Network,
    carriers: list[str] | None = None,
    horizon: int = 168,
    overlap: int = 1,
) -> None:
    """
    Fix reservoir SOC only at rolling-horizon window boundaries.

    Sets state_of_charge_set at the first and last snapshot of each RH window,
    leaving all other snapshots as NaN (unconstrained).  This guides the
    seasonal SOC trajectory while giving the optimizer freedom to choose the
    hourly dispatch path.

    Parameters
    ----------
    n : pypsa.Network
        Target network for rolling horizon (will be modified in place).
    n_msv : pypsa.Network
        Network with perfect foresight solution.
    carriers : list[str], optional
        Carriers to fix. Defaults to ["hydro-reservoir"].
    horizon : int
        Number of snapshots per RH window. Default 168 (one week at 1H).
    overlap : int
        Number of overlapping snapshots between consecutive windows. Default 1.
    """
    if carriers is None:
        carriers = ["hydro-reservoir"]

    if n_msv.storage_units_t.state_of_charge.empty:
        logger.warning(
            "PF network has no state_of_charge data, skipping SOC boundary fix"
        )
        return

    sus_i = n.storage_units[n.storage_units.carrier.isin(carriers)].index
    common = sus_i.intersection(n_msv.storage_units_t.state_of_charge.columns)

    if len(common) == 0:
        logger.info(
            f"No StorageUnits with carriers {carriers} found, skipping SOC boundary fix"
        )
        return

    pf_soc = n_msv.storage_units_t.state_of_charge[common]

    # Resample if snapshots differ
    if not n.snapshots.equals(n_msv.snapshots):
        pf_soc = resample_msv_to_target(pf_soc, n.snapshots, method="ffill")

    # Compute window boundary indices (same logic as optimize_with_rolling_horizon)
    n_sns = len(n.snapshots)
    step = horizon - overlap
    boundary_idx = set()
    for start in range(0, n_sns, step):
        end = min(n_sns - 1, start + horizon - 1)
        boundary_idx.add(start)
        boundary_idx.add(end)
    boundary_snapshots = n.snapshots[sorted(boundary_idx)]

    # Build sparse SOC set: NaN everywhere, PF values at boundaries only
    soc_sparse = pd.DataFrame(
        np.nan,
        index=n.snapshots,
        columns=common,
    )
    soc_sparse.loc[boundary_snapshots, common] = pf_soc.loc[boundary_snapshots, common]

    n.storage_units_t.state_of_charge_set[common] = soc_sparse

    stats = summarize_counts(n.storage_units.loc[common, "carrier"])
    logger.info(
        f"Fixed state_of_charge_set at {len(boundary_snapshots)} boundary snapshots "
        f"(out of {n_sns}) for {len(common)} StorageUnits:\n{stats}"
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
    accumulator_carriers = snakemake.params.get("accumulator_carriers", [])

    # Disable cyclicity (cyclic_carriers remain cyclic)
    disable_store_cyclicity(
        n,
        cyclic_carriers=cyclic_carriers,
        seasonal_carriers=seasonal_carriers,
        accumulator_carriers=accumulator_carriers,
    )

    disable_global_constraints(n)

    # disable volume limits
    disable_volume_limits(n)

    # Load MSV network and apply to configured carriers
    msv_network_path = snakemake.input.get("network_msv", None)
    if msv_network_path:
        n_msv = pypsa.Network(msv_network_path)
        resample_method = snakemake.params.get("msv_resample_method", "ffill")

        apply_msv_to_network(
            n, n_msv, seasonal_carriers, accumulator_carriers, resample_method
        )

        # Initial state from PF for seasonal only (accumulators start empty)
        set_initial_state_from_pf(n, n_msv, seasonal_carriers)

        # Fix reservoir SOC at RH window boundaries from PF trajectory
        soc_boundary_carriers = snakemake.params.get("soc_boundary_carriers", [])
        cba_solving = (
            snakemake.config.get("cba", {}).get("solving", {}).get("options", {})
        )
        fix_reservoir_soc_at_boundaries(
            n,
            n_msv,
            carriers=soc_boundary_carriers,
            horizon=cba_solving.get("horizon", 168),
            overlap=cba_solving.get("overlap", 1),
        )

    else:
        logger.info("No MSV network provided, skipping MSV application")

    # Save simplified network
    n.export_to_netcdf(snakemake.output.network)
