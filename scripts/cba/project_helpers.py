# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Helper functions for CBA project processing.
"""

import logging

logger = logging.getLogger(__name__)


def get_project_status_for_horizon(
    project, horizon: int, planning_horizons: list[int]
) -> bool:
    """
    Check if a project should be in the reference for a given planning horizon.

    Reference year to horizon mapping:
    - reference2030 → horizon 2030
    - reference2035 → closest horizon >= 2035 (typically 2040)

    Cumulative logic:
    - Horizon 2030: Include if in_reference2030=True
    - Horizon 2040: Include if in_reference2030=True OR in_reference2035=True

    Args:
        project: Project row (pd.Series) with in_referenceYYYY columns
        horizon: Target planning horizon
        planning_horizons: Available planning horizons from config

    Returns:
        True if project should be in reference for this horizon
    """
    planning_horizons = sorted(planning_horizons)

    # Check all reference year columns
    for col in project.index:
        if col.startswith("in_reference") and bool(project[col]):
            # Extract reference year (e.g., 2030 from "in_reference2030")
            ref_year = int(col.replace("in_reference", ""))

            # Map reference year to planning horizon
            if ref_year in planning_horizons:
                mapped_horizon = ref_year
            else:
                # Find closest horizon >= reference_year
                future = [h for h in planning_horizons if h >= ref_year]
                if not future:
                    continue
                mapped_horizon = min(future)

            # Include if mapped horizon <= target horizon
            if mapped_horizon <= horizon:
                return True

    return False


def log_horizon_mapping(planning_horizons: list[int]):
    """Log the reference year to planning horizon mapping for debugging."""
    logger.info("\nReference year → Planning horizon mapping:")
    for ref_year in [2030, 2035, 2040]:
        planning_horizons_sorted = sorted(planning_horizons)
        if ref_year in planning_horizons_sorted:
            logger.info(f"  reference{ref_year} → horizon {ref_year}")
        else:
            future = [h for h in planning_horizons_sorted if h >= ref_year]
            if future:
                logger.info(f"  reference{ref_year} → horizon {min(future)}")
            else:
                logger.info(f"  reference{ref_year} → (no suitable horizon)")
