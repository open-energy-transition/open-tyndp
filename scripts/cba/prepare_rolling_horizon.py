# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Prepare network for cost-benefit analysis rolling horizon dispatch.

Modifications applied:
- cyclic_carriers: keep cyclicity (short-term, e.g., batteries)
- All other stores/storage units: disable cyclicity, apply marginal storage
  value, and set initial state of charge from perfect foresight for those
  that were cyclic in the full-year optimisation
- Remove global constraints not needed for cost-benefit analysis
"""

import logging

import numpy as np
import pandas as pd
import pypsa
from numpy import inf, isfinite

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def disable_global_constraints(n: pypsa.Network):
    """
    Remove global constraints not needed for rolling horizon dispatch.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify
    """
    if "co2_sequestration_limit" in n.global_constraints.index:
        n.remove("GlobalConstraint", "co2_sequestration_limit")

    if "unsustainable biomass limit" in n.global_constraints.index:
        n.remove("GlobalConstraint", "unsustainable biomass limit")


def disable_store_cyclicity(
    n: pypsa.Network,
    cyclic_carriers: list[str] | None = None,
):
    """
    Enforce cyclic state-of-charge only for cyclic_carriers within each
    rolling horizon window. All other stores and storage units are made
    non-cyclic for long-term seasonal storage.

    Parameters
    ----------
    n : pypsa.Network
        Network to modify in place.
    cyclic_carriers : list[str], optional
        Carriers that remain cyclic. Defaults to empty list.
    """
    if cyclic_carriers is None:
        cyclic_carriers = []

    # Disable cyclicity for stores (except cyclic_carriers)
    has_e_cyclic = n.stores["e_cyclic"]
    is_cyclic_carrier = n.stores["carrier"].isin(cyclic_carriers)
    to_disable = has_e_cyclic & ~is_cyclic_carrier

    if to_disable.any():
        n.stores.loc[to_disable, "e_cyclic"] = False

    # Disable cyclicity for storage units (except cyclic_carriers)
    has_cyclic_soc = n.storage_units["cyclic_state_of_charge"]
    is_cyclic_carrier_su = n.storage_units["carrier"].isin(cyclic_carriers)
    to_disable_su = has_cyclic_soc & ~is_cyclic_carrier_su

    if to_disable_su.any():
        n.storage_units.loc[to_disable_su, "cyclic_state_of_charge"] = False


def resample_msv_to_target(
    msv: pd.DataFrame,
    target_snapshots: pd.DatetimeIndex,
    method: str = "ffill",
) -> pd.DataFrame:
    """
    Resample marginal storage value from extraction resolution to target resolution.

    Parameters
    ----------
    msv : pd.DataFrame
        MSV data from extraction (e.g., 24H resolution)
    target_snapshots : pd.DatetimeIndex
        Target snapshots (e.g., 3H resolution)
    method : str, optional
        Resampling method:
        - "ffill": Forward fill - each MSV value applies until the next one
        - "interpolate": Linear interpolation between marginal storage values
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
            c.static.loc[has_e_sum_min, "e_sum_min"] = -inf
            c.static.loc[has_e_sum_min, "e_sum_max"] = inf
            c.dynamic.p_set.loc[:, c.static[has_e_sum_min].index] = c.dynamic.p.loc[
                :, c.static[has_e_sum_min].index
            ]


def apply_msv_to_network(
    n: pypsa.Network,
    n_msv: pypsa.Network,
    cyclic_carriers: list[str],
    resample_method: str = "ffill",
) -> None:
    """
    Apply marginal storage value as marginal_cost to all Store and
    StorageUnit components whose carrier is NOT in cyclic_carriers.

    Parameters
    ----------
    n : pypsa.Network
        Target network (will be modified in place).
    n_msv : pypsa.Network
        Network with perfect foresight solution containing mu_energy_balance.
    cyclic_carriers : list[str]
        Carriers that remain cyclic (do not receive MSV).
    resample_method : str, optional
        Method for resampling MSV to target resolution. Default "ffill".
    """
    all_msv_carriers = [
        c
        for c in list(n.stores.carrier.unique())
        + list(n.storage_units.carrier.unique())
        if c not in cyclic_carriers
    ]

    if not all_msv_carriers:
        return

    if n_msv.stores_t.mu_energy_balance.empty:
        logger.warning(
            "Network has no mu_energy_balance. "
            "Ensure extraction was solved with assign_all_duals=True."
        )
        return

    # Check if resampling is needed
    needs_resample = not n.snapshots.equals(n_msv.snapshots)

    for c in ["Store", "StorageUnit"]:
        # Filter to non-cyclic carriers
        s_i = n_msv.c[c].static[n_msv.c[c].static.carrier.isin(all_msv_carriers)].index
        if s_i.empty:
            continue
        # Get shadow prices
        msv = n_msv.c[c].dynamic["mu_energy_balance"][s_i]
        # Resample if needed
        if needs_resample:
            msv = resample_msv_to_target(msv, n.snapshots, method=resample_method)
        # Set shadow prices as marginal cost
        n.c[c].dynamic["marginal_cost"].loc[:, s_i] = msv


def set_initial_state_from_pf(
    n: pypsa.Network,
    n_msv: pypsa.Network,
    cyclic_carriers: list[str],
) -> None:
    """
    Set initial storage state from perfect foresight solution.

    For stores and storage units that are NOT in cyclic_carriers but WERE
    cyclic in the full-year optimisation, set their initial state of charge
    to the value at the last snapshot of the perfect foresight solution.

    Must be called BEFORE disable_store_cyclicity so that the original
    cyclicity flags are still available.

    Parameters
    ----------
    n : pypsa.Network
        Target network (will be modified in place).
    n_msv : pypsa.Network
        Network with perfect foresight solution.
    cyclic_carriers : list[str]
        Carriers that remain cyclic in rolling horizon windows.
    """
    # Stores: set e_initial for non-cyclic carriers that were cyclic
    if not n_msv.stores_t.e.empty:
        pf_e_final = n_msv.stores_t.e.iloc[-1]
        was_cyclic = n.stores["e_cyclic"]
        is_cyclic_carrier = n.stores["carrier"].isin(cyclic_carriers)
        needs_initial = was_cyclic & ~is_cyclic_carrier
        candidates = n.stores[needs_initial].index.intersection(pf_e_final.index)
        if len(candidates) > 0:
            n.stores.loc[candidates, "e_initial"] = pf_e_final.loc[candidates]

    # Storage units: set state_of_charge_initial for non-cyclic carriers
    if not n_msv.storage_units_t.state_of_charge.empty:
        pf_soc_final = n_msv.storage_units_t.state_of_charge.iloc[-1]
        was_cyclic_su = n.storage_units["cyclic_state_of_charge"]
        is_cyclic_carrier_su = n.storage_units["carrier"].isin(cyclic_carriers)
        needs_initial_su = was_cyclic_su & ~is_cyclic_carrier_su
        candidates_su = n.storage_units[needs_initial_su].index.intersection(
            pf_soc_final.index
        )
        if len(candidates_su) > 0:
            n.storage_units.loc[candidates_su, "state_of_charge_initial"] = (
                pf_soc_final.loc[candidates_su]
            )


def fix_reservoir_soc_at_boundaries(
    n: pypsa.Network,
    n_msv: pypsa.Network,
    carriers: list[str] | None = None,
    horizon: int = 168,
    overlap: int = 1,
) -> None:
    """
    Fix reservoir state of charge at rolling horizon window boundaries.

    Sets state_of_charge_set at the first and last snapshot of each window,
    leaving all other snapshots unconstrained. This guides the seasonal
    trajectory while giving the optimizer freedom for hourly dispatch.

    Parameters
    ----------
    n : pypsa.Network
        Target network for rolling horizon (will be modified in place).
    n_msv : pypsa.Network
        Network with perfect foresight solution.
    carriers : list[str], optional
        Carriers to fix. Defaults to ["hydro-reservoir"].
    horizon : int
        Number of snapshots per rolling horizon window. Default 168 (one week at 1H).
    overlap : int
        Number of overlapping snapshots between consecutive windows. Default 1.
    """
    if carriers is None:
        carriers = ["hydro-reservoir"]

    if n_msv.storage_units_t.state_of_charge.empty:
        logger.warning(
            "Perfect foresight network has no state_of_charge data, "
            "skipping state of charge boundary fix"
        )
        return

    sus_i = n.storage_units[n.storage_units.carrier.isin(carriers)].index
    common = sus_i.intersection(n_msv.storage_units_t.state_of_charge.columns)

    if len(common) == 0:
        logger.debug(
            f"No StorageUnits with carriers {carriers} found, "
            "skipping state of charge boundary fix"
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

    # Build sparse state of charge set: NaN everywhere, perfect foresight values at boundaries only
    soc_sparse = pd.DataFrame(
        np.nan,
        index=n.snapshots,
        columns=common,
    )
    soc_sparse.loc[boundary_snapshots, common] = pf_soc.loc[boundary_snapshots, common]

    n.storage_units_t.state_of_charge_set[common] = soc_sparse


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

    # Load the reference network (already has fixed capacities and hurdle costs)
    n = pypsa.Network(snakemake.input.network)

    # Get storage carrier settings from config
    cyclic_carriers = snakemake.params.get("cyclic_carriers", [])

    # Load perfect foresight network
    n_msv = pypsa.Network(snakemake.input.network_msv)
    resample_method = snakemake.params.get("msv_resample_method", "ffill")

    # Set initial state from perfect foresight BEFORE disabling cyclicity
    # (needs original e_cyclic / cyclic_state_of_charge flags)
    set_initial_state_from_pf(n, n_msv, cyclic_carriers)

    # Disable cyclicity (cyclic_carriers remain cyclic)
    disable_store_cyclicity(n, cyclic_carriers=cyclic_carriers)

    # Remove global constraints
    disable_global_constraints(n)

    # Disable volume limits
    disable_volume_limits(n)

    # Apply marginal storage value to all non-cyclic carriers
    apply_msv_to_network(n, n_msv, cyclic_carriers, resample_method)

    # Fix reservoir state of charge at window boundaries from perfect foresight
    soc_boundary_carriers = snakemake.params.get("soc_boundary_carriers", [])
    cba_solving = snakemake.config.get("cba", {}).get("solving", {}).get("options", {})
    fix_reservoir_soc_at_boundaries(
        n,
        n_msv,
        carriers=soc_boundary_carriers,
        horizon=cba_solving.get("horizon", 168),
        overlap=cba_solving.get("overlap", 1),
    )

    # Save prepared network
    n.export_to_netcdf(snakemake.output.network)
