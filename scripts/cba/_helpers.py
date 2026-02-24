# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

import re

import pandas as pd


def annuity(lifetime: float, discount_rate: float) -> float:
    """Annuity factor for given lifetime (years) and discount rate."""
    if discount_rate == 0:
        return 1.0 / lifetime
    return discount_rate / (1.0 - 1.0 / (1.0 + discount_rate) ** lifetime)


def get_link_attrs(project: pd.Series, p_nom: float, annuity_factor: float) -> dict:
    """
    Return length, underwater_fraction, and capital_cost for a new DC link.

    Parameters
    ----------
    project : pd.Series
        Row from transmission_projects with columns length_km,
        underwater_fraction, capex_meur.
    p_nom : float
        Nominal capacity (MW) of this link direction.
    annuity_factor : float
        Annuity factor for annualizing the total CAPEX.
    """
    length = float(project.get("length_km", 0))
    uf = float(project.get("underwater_fraction", 0))
    capex = float(project.get("capex_meur", 0))
    length = 0.0 if pd.isna(length) else length
    uf = 0.0 if pd.isna(uf) else uf
    capex = 0.0 if pd.isna(capex) else capex
    capital_cost = capex * 1e6 * annuity_factor / p_nom if p_nom > 0 and capex > 0 else 0.0
    return dict(length=length, underwater_fraction=uf, capital_cost=capital_cost)


def filter_projects_by_specs(
    project_list: list[str], spec_list: list[str] | str | None
) -> list[str]:
    """
    Filter projects based on specifications with inclusions and exclusions.

    Supports:
    - Single projects: 't1', 's4'
    - Ranges: 't20-t25', 's4-s6'
    - Exclusions: '-t22', '-s5'
    - Exclusion ranges: '-t22-t25'

    The function operates in two modes:
    - Inclusion mode (default): Start with empty set, add specified projects
    - Removal mode: Start with all projects, remove specified ones (when first spec starts with '-')

    Parameters
    ----------
    project_list : list[str]
        List of all available project names to filter from
    spec_list : list[str], str, or None
        List of specifications, a single specification string, or None to return all projects

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
        else:
            # Single project
            op(spec)

    if not removals:
        return [p for p in project_list if p in projects]
    else:
        return [p for p in project_list if p not in projects]
