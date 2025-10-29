# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Calculate CBA indicators by comparing reference and project scenarios.

This script computes the B1 indicator (Total System Cost difference) and other
CBA metrics by analyzing the solved networks for reference and project cases.

PINT (Put In at a Time):
    - Reference: Network WITHOUT any projects
    - Project: Network WITH the specific project added
    - B1 = Cost(reference) - Cost(with project)

TOOT (Take Out One at a Time):
    - Reference: Network WITH all projects (current plan)
    - Project: Network WITHOUT the specific project (removed)
    - B1 = Cost(without project) - Cost(reference)


References:
- CBA guidelines: https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/news/2024/entso-e_4th_CBA_Guideline_240409.pdf
    - section 3.2.2: TOOT and PINT, page 23-24
- CBA implementation guidelines: https://eepublicdownloads.blob.core.windows.net/public-cdn-container/tyndp-documents/TYNDP2024/foropinion/CBA_Implementation_Guidelines.pdf
    - section 5.1: B1 - SEW, page 58-59

"""

import logging

import pandas as pd
import pypsa

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def calculate_total_system_cost(n):
    """
    Calculate total annualized system cost using PyPSA built-in statistics.

    This implementation handles:
    - Annualized capital costs
    - Time-aggregated operational costs
    - All component types (generators, links, storage, etc.)

    Args:
        n: PyPSA Network (must be solved)

    Returns:
        float: Total system cost in currency units (Euros)
    """
    if not n.is_solved:
        raise ValueError("Network must be solved before calculating costs")

    # Use PyPSA's built-in statistics methods
    capex = n.statistics.capex().sum()
    opex = n.statistics.opex(aggregate_time="sum").sum()
    total = capex + opex
    return {
        "total": total,
        "capex": capex,
        "opex": opex,
    }


def calculate_b1_indicator(n_reference, n_project, method="pint"):
    """
    Calculate B1 indicator: change in total system cost.

    The interpretation depends on the method:
    - PINT: positive B1 means beneficial (project reduces costs)
    - TOOT: positive B1 means beneficial (removing project increases costs)

    Args:
        n_reference: Reference network
        n_project: Project network
        method: Either "pint" or "toot" (case-insensitive)

    Returns:
        dict: Dictionary with B1 and component costs
    """
    method = method.lower()
    if method not in ["pint", "toot"]:
        raise ValueError(f"Method must be 'pint' or 'toot', got: {method}")

    # Calculate costs for both scenarios
    cost_reference = calculate_total_system_cost(n_reference)
    cost_project = calculate_total_system_cost(n_project)

    if method == "pint":
        # PINT: positive B1 means beneficial (project reduces costs)
        # Reference is without project
        # Project is with project
        b1 = cost_reference["total"] - cost_project["total"]
    else:  # toot
        # TOOT: positive B1 means beneficial (removing project increases costs)
        # Reference is with all projects
        # Project is without project
        b1 = cost_project["total"] - cost_reference["total"]

    is_beneficial = b1 > 0

    if method == "pint":
        # Speak in terms of sew.
        if is_beneficial:
            interpretation = "The project reduces costs compared to the reference scenario without the project."
        else:
            interpretation = "The project increases costs compared to the reference scenario without the project."
    else:  # toot
        if is_beneficial:
            interpretation = "The project is beneficial as removing it increases costs compared to the reference scenario with all projects."
        else:
            interpretation = "The project increases costs as removing it decreases costs compared to the reference scenario with all projects."

    results = {}
    results["B1_total_system_cost_change"] = b1  # in Euros. Positive is beneficial.
    results["is_beneficial"] = is_beneficial
    results["interpretation"] = interpretation
    results["cost_reference"] = cost_reference["total"]
    results["capex_reference"] = cost_reference["capex"]
    results["opex_reference"] = cost_reference["opex"]
    results["cost_project"] = cost_project["total"]
    results["capex_project"] = cost_project["capex"]
    results["opex_project"] = cost_project["opex"]

    if method == "pint":
        results["capex_change"] = cost_reference["capex"] - cost_project["capex"]
        results["opex_change"] = cost_reference["opex"] - cost_project["opex"]
    else:  # toot
        results["capex_change"] = cost_project["capex"] - cost_reference["capex"]
        results["opex_change"] = cost_project["opex"] - cost_reference["opex"]

    return results


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("make_indicators")

    # Load both networks
    n_reference = pypsa.Network(snakemake.input.reference)
    n_project = pypsa.Network(snakemake.input.project)

    # Validate networks are solved
    if not n_reference.is_solved:
        raise ValueError("Reference network is not solved")
    if not n_project.is_solved:
        raise ValueError("Project network is not solved")

    # Detect method from wildcards (toot or pint)
    method = snakemake.wildcards.cba_method.lower()

    # Calculate indicators
    indicators = calculate_b1_indicator(n_reference, n_project, method=method)

    # Add project metadata
    indicators["project_id"] = snakemake.wildcards.cba_project
    indicators["planning_horizon"] = snakemake.wildcards.planning_horizons

    # Convert to DataFrame and save
    df = pd.DataFrame([indicators])
    df.to_csv(snakemake.output.indicators, index=False)

    # Print detailed results

    logger.info(f"\n{'=' * 80}")
    logger.info("B1 INDICATOR RESULTS")
    logger.info(f"{'=' * 80}")
    logger.info(f"Project:         {snakemake.wildcards.cba_project}")
    logger.info(f"Method:          {method.upper()}")
    logger.info(f"Planning Horizon: {snakemake.wildcards.planning_horizons}")

    logger.info("\nREFERENCE SCENARIO:")
    logger.info(f"  Total Cost:  {indicators['cost_reference'] / 1e9:>12.2f} B€")
    logger.info(f"    CAPEX:     {indicators['capex_reference'] / 1e9:>12.2f} B€")
    logger.info(f"    OPEX:      {indicators['opex_reference'] / 1e9:>12.2f} B€")

    logger.info("\nPROJECT SCENARIO:")
    logger.info(f"  Total Cost:  {indicators['cost_project'] / 1e9:>12.2f} B€")
    logger.info(f"    CAPEX:     {indicators['capex_project'] / 1e9:>12.2f} B€")
    logger.info(f"    OPEX:      {indicators['opex_project'] / 1e9:>12.2f} B€")

    logger.info("\nCOST CHANGES:")
    logger.info(f"  CAPEX Δ:    {indicators['capex_change'] / 1e9:>12.2f} B€")
    logger.info(f"  OPEX Δ:     {indicators['opex_change'] / 1e9:>12.2f} B€")

    logger.info(
        f"\nB1 INDICATOR: {indicators['B1_total_system_cost_change'] / 1e9:>12.2f} B€"
    )
    # Print interpretation
    if indicators["is_beneficial"]:
        logger.info("PROJECT IS BENEFICIAL")
    else:
        logger.info("PROJECT INCREASES COSTS")

    logger.info(f"\nInterpretation: {indicators['interpretation']}")
    logger.info(f"{'=' * 80}\n")
