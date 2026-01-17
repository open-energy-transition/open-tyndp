# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
#

import pandas as pd

from scripts.cba._helpers import (
    filter_projects_by_method,
    filter_projects_by_specs,
    load_method_assignment,
)
from scripts._helpers import fill_wildcards


wildcard_constraints:
    cba_project=r"(s|t)\d+",


rule retrieve_tyndp_cba_projects:
    params:
        # TODO Integrate into Zenodo tyndp data bundle
        url="https://storage.googleapis.com/open-tyndp-data-store/CBA_projects.zip",
        source="CBA project explorer",
    input:
        "data/tyndp_2024_bundle",
    output:
        dir=directory("data/tyndp_2024_bundle/cba_projects"),
    log:
        "logs/retrieve_tyndp_cba_projects",
    retries: 2
    script:
        "../scripts/sb/retrieve_additional_tyndp_data.py"


# read in transmission and storage projects from excel sheets
#
def input_clustered_network(w):
    scenario = config_provider("scenario")(w)
    (clusters,) = scenario["clusters"]
    return fill_wildcards(rules.cluster_network.output.network, clusters=clusters)


checkpoint clean_projects:
    params:
        method_assignment=config_provider("cba", "method_assignment"),
    input:
        dir="data/tyndp_2024_bundle/cba_projects",
        network=input_clustered_network,
    output:
        transmission_projects=resources("cba/transmission_projects.csv"),
        storage_projects=resources("cba/storage_projects.csv"),
        method_assignment=resources("cba/method_assignment.csv"),
    script:
        "../scripts/cba/clean_projects.py"


def input_sb_network(w):
    scenario = config_provider("scenario")(w)
    expanded_wildcards = {
        "clusters": scenario["clusters"],
        "opts": scenario["opts"],
        "sector_opts": scenario["sector_opts"],
    }
    match config_provider("foresight")(w):
        case "perfect":
            expanded_wildcards["planning_horizons"] = "all"
        case "myopic":
            pass
        case _:
            raise ValueError('config["foresight"] must be one of "perfect" or "myopic"')

    return fill_wildcards(
        RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
        **expanded_wildcards,
    )


# extract a planning horizon from the SB optimized network and apply the simplifications
# necessary to get to the general CBA reference network
rule simplify_sb_network:
    params:
        hurdle_costs=config_provider("cba", "hurdle_costs"),
    input:
        network=input_sb_network,
    output:
        network=resources("cba/networks/simple_{planning_horizons}.nc"),
    script:
        "../scripts/cba/simplify_sb_network.py"


# build the unified CBA reference network
# Adds TOOT projects (not in base grid) to create reference for both TOOT and PINT
rule prepare_cba_reference:
    params:
        hurdle_costs=config_provider("cba", "hurdle_costs"),
    input:
        network=rules.simplify_sb_network.output.network,
        transmission_projects=rules.clean_projects.output.transmission_projects,
        storage_projects=rules.clean_projects.output.storage_projects,
        method_assignment=rules.clean_projects.output.method_assignment,
    output:
        network=resources("cba/networks/reference_{planning_horizons}.nc"),
    script:
        "../scripts/cba/prepare_cba_reference.py"


# remove the single project {cba_project} from the CBA reference network (TOOT)
rule prepare_toot_project:
    params:
        hurdle_costs=config_provider("cba", "hurdle_costs"),
    input:
        network=rules.prepare_cba_reference.output.network,
        transmission_projects=rules.clean_projects.output.transmission_projects,
        storage_projects=rules.clean_projects.output.storage_projects,
    output:
        network=resources(
            "cba/toot/networks/project_{cba_project}_{planning_horizons}.nc"
        ),
    script:
        "../scripts/cba/prepare_toot_project.py"


# add the single project {cba_project} to the CBA reference network (PINT)
rule prepare_pint_project:
    params:
        hurdle_costs=config_provider("cba", "hurdle_costs"),
    input:
        network=rules.prepare_cba_reference.output.network,
        transmission_projects=rules.clean_projects.output.transmission_projects,
        storage_projects=rules.clean_projects.output.storage_projects,
    output:
        network=resources(
            "cba/pint/networks/project_{cba_project}_{planning_horizons}.nc"
        ),
    script:
        "../scripts/cba/prepare_pint_project.py"


# solve any of the prepared networks, ie a reference or a project network
# should reuse/import functions from solve_network.py
rule solve_cba_network:
    params:
        solving=config_provider("solving"),
        cba_solving=config_provider("cba", "solving"),
        foresight=config_provider("foresight"),
        time_resolution=config_provider("clustering", "temporal", "resolution_sector"),
        custom_extra_functionality=None,
    input:
        network=resources("cba/{cba_method}/networks/{name}_{planning_horizons}.nc"),
    output:
        network=RESULTS + "cba/{cba_method}/networks/{name}_{planning_horizons}.nc",
    log:
        solver=logs(
            "cba/solve_cba_network/{cba_method}_{name}_{planning_horizons}_solver.log"
        ),
        memory=logs(
            "cba/solve_cba_network/{cba_method}_{name}_{planning_horizons}_memory.log"
        ),
        python=logs(
            "cba/solve_cba_network/{cba_method}_{name}_{planning_horizons}_python.log"
        ),
    threads: 1
    script:
        "../scripts/cba/solve_cba_network.py"


# solve the unified reference network (shared by TOOT and PINT)
rule solve_cba_reference:
    params:
        solving=config_provider("solving"),
        cba_solving=config_provider("cba", "solving"),
        foresight=config_provider("foresight"),
        time_resolution=config_provider("clustering", "temporal", "resolution_sector"),
        custom_extra_functionality=None,
    input:
        network=resources("cba/networks/reference_{planning_horizons}.nc"),
    output:
        network=RESULTS + "cba/networks/reference_{planning_horizons}.nc",
    log:
        solver=logs("cba/solve_cba_network/reference_{planning_horizons}_solver.log"),
        memory=logs("cba/solve_cba_network/reference_{planning_horizons}_memory.log"),
        python=logs("cba/solve_cba_network/reference_{planning_horizons}_python.log"),
    threads: 1
    script:
        "../scripts/cba/solve_cba_network.py"


# compute all metrics for a single pint or toot project comparing reference and project solution
rule make_indicators:
    input:
        reference=RESULTS + "cba/networks/reference_{planning_horizons}.nc",
        project=RESULTS
        + "cba/{cba_method}/networks/project_{cba_project}_{planning_horizons}.nc",
    output:
        indicators=RESULTS
        + "cba/{cba_method}/project_{cba_project}_{planning_horizons}.csv",
    script:
        "../scripts/cba/make_indicators.py"


def input_indicators(w):
    """
    List all indicator CSV files for the given CBA method.

    Filters projects based on the method (TOOT vs PINT) using method_assignment.csv
    and user-specified project specs from config.
    """
    run = w.get("run", config_provider("run", "name")(w))
    planning_horizons = int(w.planning_horizons)
    cba_method = w.cba_method
    scenario = run  # NT, DE, etc.

    transmission_projects = pd.read_csv(
        checkpoints.clean_projects.get(run=run).output.transmission_projects
    )
    storage_projects = pd.read_csv(
        checkpoints.clean_projects.get(run=run).output.storage_projects
    )
    method_assignment = load_method_assignment(
        checkpoints.clean_projects.get(run=run).output.method_assignment
    )

    # Filter projects based on CBA method using method_assignment.csv
    transmission_projects = filter_projects_by_method(
        transmission_projects, cba_method, planning_horizons, scenario, method_assignment
    )

    cba_projects = [
        f"t{pid}" for pid in transmission_projects["project_id"].unique()
    ] + [f"s{pid}" for pid in storage_projects["project_id"].unique()]

    project_specs = config_provider("cba", "projects")(w)

    return expand(
        rules.make_indicators.output.indicators,
        cba_project=filter_projects_by_specs(cba_projects, project_specs),
        allow_missing=True,
    )


# collect the indicators for all transmission_projects into a single overview csv
rule collect_indicators:
    input:
        indicators=input_indicators,
    output:
        indicators=RESULTS + "cba/{cba_method}/indicators_{planning_horizons}.csv",
    script:
        "../scripts/cba/collect_indicators.py"


rule plot_indicators:
    params:
        plotting=config_provider("plotting"),
    input:
        indicators=rules.collect_indicators.output.indicators,
        transmission_projects=rules.clean_projects.output.transmission_projects,
    output:
        plot_dir=directory(RESULTS + "cba/{cba_method}/plots_{planning_horizons}"),
    script:
        "../scripts/cba/plot_indicators.py"


# pseudo-rule, to run enable running cba with snakemake cba --configfile config/config.tyndp.yaml
rule cba:
    input:
        lambda w: expand(
            rules.collect_indicators.output.indicators,
            cba_method=config_provider("cba", "methods")(w),
            planning_horizons=config_provider("cba", "planning_horizons")(w),
            run=config_provider("run", "name")(w),
        ),
        lambda w: expand(
            rules.plot_indicators.output.plot_dir,
            cba_method=config_provider("cba", "methods")(w),
            planning_horizons=config_provider("cba", "planning_horizons")(w),
            run=config_provider("run", "name")(w),
        ),
