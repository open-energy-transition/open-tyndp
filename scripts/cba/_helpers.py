# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

import fnmatch
import re

import pandas as pd

from scripts.add_electricity import calculate_annuity


def get_storage_attrs(project: pd.Series, discount_rate: float) -> dict:
    """
    Return sized bus/link/store attributes for a new CBA storage project.

    All costs and capacities are taken directly from the CBA storage projects
    Excel export (via ``project``), rather than from the technology cost
    tables, since project-level CAPEX and OPEX are reported explicitly.
    CAPEX is annualized using the project's operational lifetime (already
    cleaned up in clean_projects.py, i.e. missing or zero values replaced
    with default_lifetime) and ``discount_rate``; the already-annual OPEX is
    added on top (matching how `capital_cost` is built from investment + FOM
    in process_cost_data.py). The combined annualized cost is booked entirely
    on the Store (in EUR/MWh/yr) so it is only counted once; both links get
    zero capital cost.

    Parameters
    ----------
    project : pd.Series
        Row from storage_projects with columns p_nom_discharge, p_nom_charge,
        e_nom_gwh, roundtrip_efficiency, capex_meur, opex_meur_per_year, and
        operational lifetime in years under "lifetime_years".
    discount_rate : float
        Discount rate used to annualize the project's CAPEX.

    Returns
    -------
    dict
        Dictionary with keys p_nom_discharge, p_nom_charge, e_nom,
        efficiency, and capital_cost_per_mwh.
    """
    annuity = calculate_annuity(project["lifetime_years"], discount_rate)
    annualized_capex = project["capex_meur"] * 1e6 * annuity
    annual_opex = project["opex_meur_per_year"] * 1e6
    annualized_cost = annualized_capex + annual_opex

    p_nom_discharge = float(project["p_nom_discharge"])
    p_nom_charge = float(project["p_nom_charge"])
    e_nom = float(project["e_nom_gwh"]) * 1e3  # GWh -> MWh

    capital_cost_per_mwh = annualized_cost / e_nom if e_nom > 0 else 0.0

    efficiency = project["roundtrip_efficiency"] ** 0.5

    return dict(
        p_nom_discharge=p_nom_discharge,
        p_nom_charge=p_nom_charge,
        e_nom=e_nom,
        efficiency=efficiency,
        capital_cost_per_mwh=capital_cost_per_mwh,
    )


def get_link_attrs(project: pd.Series, costs: pd.DataFrame) -> dict:
    """
    Return length, underwater_fraction, and capital_cost for a new DC link.

    The capital_cost is computed using the same per-km formula as
    `add_electricity.py` to ensure consistency with existing network
    links:

    `capital_cost = length * ((1 - uf) * overhead + uf * submarine) + inverter`

    Parameters
    ----------
    project : pd.Series
        Row from transmission_projects with columns length_km and
        underwater_fraction.
    costs : pd.DataFrame
        Technology costs table (indexed by technology name) with a
        `capital_cost` column containing annualized EUR/MW or EUR/MW/km
        values.

    Returns
    -------
    dict
        Dictionary with keys length, underwater_fraction, and capital_cost.
    """
    length = float(project.get("length_km", 0))
    uf = float(project.get("underwater_fraction", 0))
    length = 0.0 if pd.isna(length) else length
    uf = 0.0 if pd.isna(uf) else uf
    overhead = costs.at["HVDC overhead", "capital_cost"]
    submarine = costs.at["HVDC submarine", "capital_cost"]
    inverter = costs.at["HVDC inverter pair", "capital_cost"]
    capital_cost = length * ((1.0 - uf) * overhead + uf * submarine) + inverter
    return dict(length=length, underwater_fraction=uf, capital_cost=capital_cost)


def filter_projects_by_specs(
    project_list: list[str], spec_list: list[str] | str | None
) -> list[str]:
    """
    Filter projects based on specifications with inclusions and exclusions.

    Supports:
    - Single projects: 't1', 's4'
    - Ranges: 't20-t25', 's4-s6'
    - Globs: 's*', 't1*'
    - Exclusions: '-t22', '-s5'
    - Exclusion ranges: '-t22-t25'
    - Exclusion globs: '-s*'

    The function operates in two modes:
    - Inclusion mode (default): Start with empty set, add specified projects
    - Removal mode: Start with all projects, remove specified ones (when first spec starts with '-')

    Parameters
    ----------
    project_list : list[str]
        List of all available project names to filter from.
    spec_list : list[str], str, or None
        List of specifications, a single specification string, or None to return all projects.

    Returns
    -------
    list[str]
        Filtered list of projects, preserving order from project_list

    Examples
    --------
    >>> filter_projects_by_specs(['t20', 't21', 't22', 't23'], ['t20-t22'])
    ['t20', 't21', 't22']

    >>> filter_projects_by_specs(['t20', 't21', 't22', 't23'], 't20-t22')
    ['t20', 't21', 't22']

    >>> filter_projects_by_specs(['t20', 't21', 't22', 't23'], ['t20-t23', '-t22'])
    ['t20', 't21', 't23']

    >>> filter_projects_by_specs(['t1', 't2', 't3', 't4'], ['-t2', '-t3'])
    ['t1', 't4']

    >>> filter_projects_by_specs(['t1', 's1', 's2'], ['s*'])
    ['s1', 's2']
    """

    if not spec_list:
        return project_list

    if isinstance(spec_list, str):
        spec_list = [spec_list]

    projects = set()
    range_pattern = re.compile(r"^([a-z])(\d+)-\1(\d+)$")

    removals = spec_list[0].startswith("-")

    for spec in spec_list:
        # Check if this is an exclusion
        if spec.startswith("-"):
            spec = spec[1:]
            op = projects.discard if not removals else projects.add
        else:
            op = projects.add if not removals else projects.discard

        # Try to match range pattern
        match = range_pattern.match(spec)
        if match:
            prefix = match.group(1)
            start = int(match.group(2))
            end = int(match.group(3))

            for i in range(start, end + 1):
                op(f"{prefix}{i}")
        elif any(c in spec for c in "*?["):
            # Glob pattern
            for p in fnmatch.filter(project_list, spec):
                op(p)
        else:
            # Single project
            op(spec)

    if not removals:
        filtered_list = [p for p in project_list if p in projects]
    else:
        filtered_list = [p for p in project_list if p not in projects]

    if not filtered_list:
        raise ValueError(f"Project specification {spec_list} selects no projects.")

    return filtered_list
