# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
#

import os
import pandas as pd

from scripts.cba._helpers import filter_projects_by_specs
from scripts._helpers import fill_wildcards
from shutil import unpack_archive, copy2


wildcard_constraints:
    cba_project=r"(s|t)\d+",


if (CBA_PROJECTS_DATASET := dataset_version("tyndp_cba_projects"))["source"] in [
    "archive"
]:

    rule retrieve_tyndp_cba_projects:
        params:
            source="CBA project explorer",
        input:
            # TODO Integrate into Zenodo tyndp data bundle
            zip_file=storage(CBA_PROJECTS_DATASET["url"]),
            dir=rules.retrieve_tyndp.output.dir,
        output:
            dir=directory(CBA_PROJECTS_DATASET["folder"]),
        log:
            "logs/retrieve_tyndp_cba_projects",
        run:
            copy2(input["zip_file"], output["dir"] + ".zip")
            unpack_archive(output["dir"] + ".zip", output["dir"])
            os.remove(output["dir"] + ".zip")


if (CBA_NON_CO2_DATASET := dataset_version("tyndp_cba_non_co2_emissions"))[
    "source"
] in ["archive"]:

    rule retrieve_tyndp_cba_non_co2_emissions:
        input:
            file=storage(CBA_NON_CO2_DATASET["url"]),
        output:
            file=resources("cba/a.3_non-co2-emissions.csv"),
        log:
            logs("retrieve_tyndp_cba_non_co2_emissions.log"),
        run:
            copy2(input["file"], output["file"])


# read in transmission and storage projects from excel sheets
#
def input_clustered_network(w):
    scenario = config_provider("scenario")(w)
    (clusters,) = scenario["clusters"]
    return fill_wildcards(rules.cluster_network.output.network, clusters=clusters)


checkpoint clean_projects:
    input:
        dir=rules.retrieve_tyndp_cba_projects.output.dir,
        network=input_clustered_network,
    output:
        transmission_projects=resources("cba/transmission_projects.csv"),
        storage_projects=resources("cba/storage_projects.csv"),
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
        tyndp_conventional_carriers=config_provider(
            "electricity", "tyndp_conventional_carriers"
        ),
        cyclic_carriers=config_provider("cba", "storage", "cyclic_carriers"),
    input:
        network=input_sb_network,
    output:
        network=resources("cba/networks/simple_{planning_horizons}.nc"),
    script:
        "../scripts/cba/simplify_sb_network.py"


# build the reference network
rule prepare_reference:
    params:
        hurdle_costs=config_provider("cba", "hurdle_costs"),
    input:
        network=rules.simplify_sb_network.output.network,
        transmission_projects=rules.clean_projects.output.transmission_projects,
        storage_projects=rules.clean_projects.output.storage_projects,
    output:
        network=resources("cba/networks/reference_{planning_horizons}.nc"),
    script:
        "../scripts/cba/prepare_reference.py"


# remove the single project {cba_project} from the toot reference network
# currently this can be either a trans123 or a stor123 project
rule prepare_toot_project:
    params:
        hurdle_costs=config_provider("cba", "hurdle_costs"),
    input:
        network=rules.prepare_reference.output.network,
        transmission_projects=rules.clean_projects.output.transmission_projects,
        storage_projects=rules.clean_projects.output.storage_projects,
    output:
        network=resources(
            "cba/toot/networks/project_{cba_project}_{planning_horizons}.nc"
        ),
    script:
        "../scripts/cba/prepare_toot_project.py"


# add the single project {cba_project} to the pint reference network
# currently this can be either a trans123 or a stor123 project
rule prepare_pint_project:
    params:
        hurdle_costs=config_provider("cba", "hurdle_costs"),
    input:
        network=rules.prepare_reference.output.network,
        transmission_projects=rules.clean_projects.output.transmission_projects,
        storage_projects=rules.clean_projects.output.storage_projects,
    output:
        network=resources(
            "cba/pint/networks/project_{cba_project}_{planning_horizons}.nc"
        ),
    script:
        "../scripts/cba/prepare_pint_project.py"


# solve the reference network which is independent of the method
rule solve_cba_reference_network:
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
        solver=logs(
            "cba/solve_cba_reference_network/reference_{planning_horizons}_solver.log"
        ),
        memory=logs(
            "cba/solve_cba_reference_network/reference_{planning_horizons}_memory.log"
        ),
        python=logs(
            "cba/solve_cba_reference_network/reference_{planning_horizons}_python.log"
        ),
    threads: 1
    script:
        "../scripts/cba/solve_cba_network.py"


# solve any of the prepared project network
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


# compute all metrics for a single pint or toot project comparing reference and project solution
rule make_indicators:
    input:
        reference=RESULTS + "cba/networks/reference_{planning_horizons}.nc",
        project=RESULTS
        + "cba/{cba_method}/networks/project_{cba_project}_{planning_horizons}.nc",
        non_co2_emissions=rules.retrieve_tyndp_cba_non_co2_emissions.output.file,
    output:
        indicators=RESULTS
        + "cba/{cba_method}/project_{cba_project}_{planning_horizons}.csv",
    script:
        "../scripts/cba/make_indicators.py"


def input_indicators(w):
    """
    List all indicators csv
    """
    run = w.get("run", config_provider("run", "name")(w))
    transmission_projects = pd.read_csv(
        checkpoints.clean_projects.get(run=run).output.transmission_projects
    )
    storage_projects = pd.read_csv(
        checkpoints.clean_projects.get(run=run).output.storage_projects
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
