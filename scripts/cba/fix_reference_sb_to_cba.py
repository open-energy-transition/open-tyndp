# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Compute the capacity corrections from the TYNDP 2026 reference grid total to the CBA reference grid
of a planning horizon, following the Implementation Guidelines Appendix Table 3.
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.sb._helpers import extract_grid_data_tyndp, format_bz_names

logger = logging.getLogger(__name__)


def read_reference_grid(grid_fn: str, horizon: int) -> pd.DataFrame:
    """
    Read the TYNDP 2026 electricity reference grid of a planning horizon per direction.

    Parameters
    ----------
    grid_fn : str
        Path to the TYNDP 2026 electricity reference grid workbook.
    horizon : int
        Planning horizon, selecting the sheet `Year_{horizon}`.

    Returns
    -------
    pd.DataFrame
        Table indexed by directed link ID (`{bus0}-{bus1}-DC`) with the 2025 NTC
        (`ntc25`), the reference grid total (`total`) and one column per tagged
        project ID ("p<id>_d" or "p<id>_in") with its capacity.
    """
    grid = pd.read_excel(grid_fn, sheet_name=f"Year_{horizon}")
    grid["Border"] = format_bz_names(grid["Border"])
    tags = grid.columns[grid.columns.str.fullmatch(r"P\d+_d")].str[1:-2]
    value_columns = {
        "ntc25": ("TY2026_2025TY Direct NTC (MW)", "TY2026_2025TY Indirect NTC (MW)"),
        "total": ("Summary Direction 1", "Summary Direction 2"),
    } | {int(t): (f"P{t}_d", f"P{t}_in") for t in tags}
    grid = extract_grid_data_tyndp(
        grid, idx_suffix="DC", idx_separator="-", value_columns=value_columns
    )
    return grid.drop(columns=["bus0", "bus1"]).fillna(0).rename_axis("link_id")


def reconcile_reference_grid(
    grid: pd.DataFrame, projects: pd.DataFrame, t3_methods: pd.Series
) -> pd.Series:
    """
    Compute the capacity corrections from the reference grid total to the CBA reference grid.

    The CBA reference grid follows the Implementation Guidelines Appendix Table 3:

    - TOOT projects must be in the reference grid. If the reference grid total is lower than the 2025 NTC
      plus the TOOT projects on a link, the missing capacity is added.
    - PINT projects must not be in the reference grid. The capacity above the 2025 NTC and the TOOT projects
      is unassigned. In project ID order, PINT projects tagged in the reference grid and then projects not
      in Table 3 are removed from the reference grid, if their capacity is not larger than the unassigned
      capacity left on each of their links. Projects with links to offshore hubs (`*_OFF` buses) are never
      removed, as their offshore wind would be cut off from the grid.

    Parameters
    ----------
    grid : pd.DataFrame
        Reference grid per directed link, see `read_reference_grid`.
    projects : pd.DataFrame
        Transmission projects with `bus0`, `bus1`, `p_nom 0->1` and `p_nom 1->0`.
    t3_methods : pd.Series
        Method (`toot` or `pint`) per project ID from draft CBA Implementation Guidelines Appendix Table 3, for the planning horizon.

    Returns
    -------
    pd.Series
        Capacity corrections per directed link ID.
    """
    directions = [("bus0", "bus1", "p_nom 0->1"), ("bus1", "bus0", "p_nom 1->0")]
    directed = (
        pd.concat(
            projects.assign(
                link_id=projects[bus0] + "-" + projects[bus1] + "-DC",
                capacity=projects[p_nom],
            )
            for bus0, bus1, p_nom in directions
        )
        .query("capacity > 0")
        .groupby(["project_id", "link_id"])["capacity"]
        .sum()
    )

    method = directed.index.get_level_values("project_id").map(t3_methods)
    toot = directed[method == "toot"].groupby("link_id").sum()
    unassigned = (grid["total"] - grid["ntc25"]).sub(toot, fill_value=0)
    corrections = (-unassigned).clip(lower=0).reindex(toot.index)
    available = unassigned.clip(lower=0)

    tags = grid.drop(columns=["ntc25", "total"])
    tagged = tags.columns[tags.gt(0).any()]
    project_ids = directed.index.unique("project_id")
    t3_pint = t3_methods.index[t3_methods.eq("pint")].intersection(tagged)

    missing_toot = (
        t3_methods.index[t3_methods.eq("toot")]
        .intersection(tagged)
        .difference(project_ids)
    )
    if not missing_toot.empty:
        logger.warning(
            f"TOOT projects {missing_toot.tolist()} (from CBA Implementation Guide Appendix Table 3) are tagged in the reference grid but not in the cleaned transmission projects (e.g. dropped for unknown bus codes), so their capacity is treated as unassigned."
        )

    candidates = [
        *project_ids.intersection(t3_pint).sort_values(),
        *project_ids.difference(t3_methods.index).sort_values(),
    ]

    to_hub = projects[["bus0", "bus1"]].apply(lambda bus: bus.str.endswith("_OFF"))
    hub_ids = set(projects.loc[to_hub.any(axis=1), "project_id"])
    hub_projects = [i for i in candidates if i in hub_ids]
    if hub_projects:
        logger.info(
            f"Not removing projects {hub_projects} from the reference grid, as they connect offshore hubs and removing them would cut off offshore wind."
        )
    candidates = [i for i in candidates if i not in hub_projects]

    not_removed = []
    for project_id in candidates:
        capacity = directed[project_id]
        if capacity.gt(available.reindex(capacity.index, fill_value=0)).any():
            if project_id in t3_pint:
                not_removed.append(project_id)
            continue
        available = available.sub(capacity, fill_value=0)
        corrections = corrections.sub(capacity, fill_value=0)

    if not_removed:
        logger.warning(
            f"Not removing PINT projects {not_removed} from the reference grid, as their capacity is larger than the unassigned capacity on their links."
        )

    return corrections[corrections.ne(0)].rename("p_nom")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "fix_reference_sb_to_cba",
            horizon="2040",
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    horizon = int(snakemake.wildcards.horizon)
    projects = pd.read_csv(snakemake.input.transmission_projects)
    in_ref = pd.read_csv(snakemake.input.guidelines).set_index("ID")[
        f"In_ref_grid_{horizon}"
    ]
    is_toot = in_ref.str.strip().str.lower().eq("yes").groupby(level=0).any()
    t3_methods = is_toot.map({True: "toot", False: "pint"})

    grid = read_reference_grid(snakemake.input.elec_reference_grid, horizon)
    corrections = reconcile_reference_grid(grid, projects, t3_methods)
    logger.info(
        f"Reference grid corrections for {horizon} (MW):\n{corrections.to_string()}"
    )
    corrections.to_csv(snakemake.output.corrections)
