# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Assign CBA method (TOOT/PINT) per project and planning horizon.

Uses the CBA Implementation Guidelines reference grid table (Table B1)
and the cleaned transmission projects list.
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def normalize_yes_no(value: str) -> str:
    return str(value).strip().lower()


def compute_method(flag: str) -> str:
    return "TOOT" if flag == "yes" else "PINT"


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("assign_cba_project_method")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    guidelines = pd.read_csv(snakemake.input.guidelines)
    projects = pd.read_csv(snakemake.input.transmission_projects)

    guidelines = guidelines.rename(
        columns={
            "ID": "project_id",
            "Project_name": "project_name",
            "In_ref_grid_2030": "in_ref_2030",
            "In_ref_grid_2040": "in_ref_2040",
        }
    )

    for col in ["in_ref_2030", "in_ref_2040"]:
        if col in guidelines.columns:
            guidelines[col] = guidelines[col].map(normalize_yes_no)

    base = guidelines[["project_id", "project_name", "in_ref_2030", "in_ref_2040"]]
    base = base.dropna(subset=["project_id"])

    # Aggregate duplicates: any "yes" implies TOOT.
    agg = base.groupby("project_id", as_index=False).agg(
        project_name=("project_name", "first"),
        in_ref_2030=("in_ref_2030", lambda s: "yes" if (s == "yes").any() else "no"),
        in_ref_2040=("in_ref_2040", lambda s: "yes" if (s == "yes").any() else "no"),
    )

    assigned = []
    for horizon, col in [(2030, "in_ref_2030"), (2040, "in_ref_2040")]:
        rows = agg[["project_id", "project_name", "in_ref_2030", "in_ref_2040"]].copy()
        rows["planning_horizon"] = horizon
        rows["method"] = rows[col].map(compute_method)
        rows = rows.rename(
            columns={
                "in_ref_2030": "in_ref_grid_2030",
                "in_ref_2040": "in_ref_grid_2040",
            }
        )
        assigned.append(rows)

    assigned = pd.concat(assigned, ignore_index=True)

    # Join to transmission projects and warn on missing IDs.
    merged = projects.merge(
        assigned,
        on="project_id",
        how="left",
        suffixes=("", "_guideline"),
    )

    missing = merged["method"].isna()
    if missing.any():
        logger.warning(
            "Missing method assignments for %d transmission projects: %s",
            missing.sum(),
            merged.loc[missing, "project_id"].tolist(),
        )

    # TODO: extend to storage projects when storage extraction is implemented.

    output = merged.copy()
    output["method"] = output["method"].fillna("PINT")
    if "in_ref_grid_2030" in output.columns:
        output["in_ref_grid_2030"] = output["in_ref_grid_2030"].fillna("no")
    else:
        output["in_ref_grid_2030"] = "no"
    if "in_ref_grid_2040" in output.columns:
        output["in_ref_grid_2040"] = output["in_ref_grid_2040"].fillna("no")
    else:
        output["in_ref_grid_2040"] = "no"
    output.to_csv(snakemake.output.methods, index=False)
