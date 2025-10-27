"""
Calculate CBA indicators by comparing reference and project scenarios.

This script computes the B1 indicator (Total System Cost difference) and other
CBA metrics by analyzing the solved networks for reference and project cases.

TOOT (Take Out One at a Time):
    - Reference: Network WITH all projects (current plan)
    - Project: Network WITHOUT the specific project (removed)
    - B1 = Cost(without project) - Cost(with all projects)
    - Positive B1 = project is beneficial (removing it increases costs)

PINT (Put In at a Time):
    - Reference: Network WITHOUT any projects (base case)
    - Project: Network WITH the specific project added
    - B1 = Cost(with project) - Cost(without projects)
    - Negative B1 = project is beneficial (adding it reduces costs)
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


def calculate_b1_indicator(n_reference, n_project, method="pint"):
    """
    Calculate B1 indicator: change in total system cost.

    The interpretation depends on the method:

    PINT (Put In at a Time) - default:
        Reference: base case WITHOUT projects
        Project: base case WITH project added
        B1 = Cost(with project) - Cost(reference)
        - Negative B1: project reduces costs (beneficial)
        - Positive B1: project increases costs

    TOOT (Take Out One at a Time):
        Reference: current plan WITH all projects
        Project: current plan WITHOUT specific project
        B1 = Cost(without project) - Cost(reference)
        - Positive B1: project is beneficial (removing it increases costs)
        - Negative B1: project increases costs

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

    # Calculate detailed breakdown
    detail_reference = calculate_total_system_cost_detailed(n_reference)
    detail_project = calculate_total_system_cost_detailed(n_project)

    # Calculate B1 - always project minus reference
    b1 = cost_project - cost_reference

    # Determine if project is beneficial based on method
    if method == "pint":
        # PINT: negative B1 means beneficial (project reduces costs)
        is_beneficial = b1 < 0
        interpretation = "PINT: Negative B1 means beneficial (project reduces costs)"
    else:  # toot
        # TOOT: positive B1 means beneficial (removing project increases costs)
        is_beneficial = b1 > 0
        interpretation = (
            "TOOT: Positive B1 means beneficial (project prevents cost increase)"
        )

    return {
        "B1_total_system_cost_change": b1,
        "method": method.upper(),
        "is_beneficial": is_beneficial,
        "interpretation": interpretation,
        "cost_reference": cost_reference,
        "cost_project": cost_project,
        "cost_difference": abs(b1),
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
    print(f"\n{'=' * 80}")
    print("B1 INDICATOR RESULTS")
    print(f"{'=' * 80}")
    print(f"Project:         {snakemake.wildcards.cba_project}")
    print(f"Method:          {method.upper()}")
    print(f"Planning Horizon: {snakemake.wildcards.planning_horizons}")

    print("\nREFERENCE SCENARIO:")
    print(f"  Total Cost:  {indicators['cost_reference'] / 1e9:>12.2f} B€")
    print(f"    CAPEX:     {indicators['capex_reference'] / 1e9:>12.2f} B€")
    print(f"    OPEX:      {indicators['opex_reference'] / 1e9:>12.2f} B€")

    print("\nPROJECT SCENARIO:")
    print(f"  Total Cost:  {indicators['cost_project'] / 1e9:>12.2f} B€")
    print(f"    CAPEX:     {indicators['capex_project'] / 1e9:>12.2f} B€")
    print(f"    OPEX:      {indicators['opex_project'] / 1e9:>12.2f} B€")

    print("\nCOST CHANGES:")
    print(f"  CAPEX Δ:    {indicators['capex_change'] / 1e9:>12.2f} B€")
    print(f"  OPEX Δ:     {indicators['opex_change'] / 1e9:>12.2f} B€")
    print(f"  Total Δ:    {indicators['B1_total_system_cost_change'] / 1e9:>12.2f} B€")

    print(
        f"\nB1 INDICATOR: {indicators['B1_total_system_cost_change'] / 1e9:>12.2f} B€"
    )

    # Print interpretation
    if indicators["is_beneficial"]:
        print("PROJECT IS BENEFICIAL")
    else:
        print("PROJECT INCREASES COSTS")

    print(f"\nInterpretation: {indicators['interpretation']}")
    print(f"{'=' * 80}\n")
