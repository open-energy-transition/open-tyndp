"""
Calculate CBA indicators by comparing reference and project scenarios.

This script computes the B1 indicator (Total System Cost difference) and other
CBA metrics by analyzing the solved networks for reference and project cases.
"""

import pandas as pd
import pypsa


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

    return capex + opex


def calculate_total_system_cost_detailed(n):
    """
    Calculate total system cost with component breakdown.

    Returns:
        dict: Dictionary with total cost and component breakdown
    """
    if not n.is_solved:
        raise ValueError("Network must be solved before calculating costs")

    # Get CAPEX and OPEX using built-in methods
    capex_series = n.statistics.capex()
    opex_series = n.statistics.opex(aggregate_time="sum")

    # Group by component type
    capex_by_component = {}
    opex_by_component = {}

    if isinstance(capex_series.index, pd.MultiIndex):
        capex_by_component = capex_series.groupby(level=0).sum().to_dict()
    if isinstance(opex_series.index, pd.MultiIndex):
        opex_by_component = opex_series.groupby(level=0).sum().to_dict()

    total_capex = capex_series.sum()
    total_opex = opex_series.sum()
    total_cost = total_capex + total_opex

    return {
        "total_cost": total_cost,
        "total_capex": total_capex,
        "total_opex": total_opex,
        "capex_by_component": capex_by_component,
        "opex_by_component": opex_by_component,
    }


def calculate_b1_indicator(n_reference, n_project):
    """
    Calculate B1 indicator: change in total system cost.

    B1 = Cost(with project) - Cost(without project)
    - Negative B1: project reduces costs (beneficial)
    - Positive B1: project increases costs
    - Near-zero B1: project has minimal cost impact

    Args:
        n_reference: Reference network (without project)
        n_project: Project network (with project)

    Returns:
        dict: Dictionary with B1 and component costs
    """
    # Calculate costs for both scenarios
    cost_reference = calculate_total_system_cost(n_reference)
    cost_project = calculate_total_system_cost(n_project)

    # Calculate detailed breakdown
    detail_reference = calculate_total_system_cost_detailed(n_reference)
    detail_project = calculate_total_system_cost_detailed(n_project)

    # B1 indicator (positive = increased cost)
    b1 = cost_project - cost_reference

    return {
        "B1_total_system_cost_change": b1,
        "cost_reference": cost_reference,
        "cost_project": cost_project,
        "cost_savings": -b1,  # Positive means savings
        "capex_reference": detail_reference["total_capex"],
        "capex_project": detail_project["total_capex"],
        "opex_reference": detail_reference["total_opex"],
        "opex_project": detail_project["total_opex"],
        "capex_change": detail_project["total_capex"] - detail_reference["total_capex"],
        "opex_change": detail_project["total_opex"] - detail_reference["total_opex"],
    }


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

    # Calculate indicators
    indicators = calculate_b1_indicator(n_reference, n_project)

    # Add project metadata
    indicators["project_id"] = snakemake.wildcards.project

    # Convert to DataFrame and save
    df = pd.DataFrame([indicators])
    df.to_csv(snakemake.output.indicators, index=False)

    # Print detailed results
    print(f"\n{'=' * 80}")
    print(f"B1 INDICATOR RESULTS: {snakemake.wildcards.project}")
    print(f"{'=' * 80}")
    print("\nREFERENCE SCENARIO (without project):")
    print(f"  Total Cost:  {indicators['cost_reference'] / 1e9:>12.2f} B€")
    print(f"    CAPEX:     {indicators['capex_reference'] / 1e9:>12.2f} B€")
    print(f"    OPEX:      {indicators['opex_reference'] / 1e9:>12.2f} B€")

    print("\nPROJECT SCENARIO (with project):")
    print(f"  Total Cost:  {indicators['cost_project'] / 1e9:>12.2f} B€")
    print(f"    CAPEX:     {indicators['capex_project'] / 1e9:>12.2f} B€")
    print(f"    OPEX:      {indicators['opex_project'] / 1e9:>12.2f} B€")

    print("\nCOST CHANGES:")
    print(f"  CAPEX Δ:    {indicators['capex_change'] / 1e9:>12.2f} B€")
    print(f"  OPEX Δ:     {indicators['opex_change'] / 1e9:>12.2f} B€")
    print(f"  Total Δ:    {indicators['B1_total_system_cost_change'] / 1e9:>12.2f} B€")

    print(f"\nB1 INDICATOR: {indicators['B1_total_system_cost_change'] / 1e9:>12.2f} B€")

    print(f"{'=' * 80}\n")
