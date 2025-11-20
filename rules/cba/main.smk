# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
#

import pandas as pd

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
        "../../scripts/retrieve_additional_tyndp_data.py"


# read in transmission and storage projects from excel sheets
#
def input_clustered_network(w):
    scenario = config_provider("scenario")(w)
    (clusters,) = scenario["clusters"]
    return fill_wildcards(rules.cluster_network.output.network, clusters=clusters)


checkpoint clean_projects:
    input:
        dir="data/tyndp_2024_bundle/cba_projects",
        network=input_clustered_network,
    output:
        transmission_projects=resources("cba/transmission_projects.csv"),
        storage_projects=resources("cba/storage_projects.csv"),
    script:
        "../../scripts/cba/clean_projects.py"


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
    input:
        network=input_sb_network,
    output:
        network=resources("cba/networks/simple_{planning_horizons}.nc"),
    script:
        "../../scripts/cba/simplify_sb_network.py"


# build the reference network for toot with all projects included
# maybe worth to merge with pint rule, if they turn out to be very similar
rule prepare_toot_reference:
    input:
        network=rules.simplify_sb_network.output.network,
        transmission_projects=rules.clean_projects.output.transmission_projects,
        storage_projects=rules.clean_projects.output.storage_projects,
    output:
        network=resources("cba/toot/networks/reference_{planning_horizons}.nc"),
    script:
        "../../scripts/cba/prepare_toot_reference.py"


# build the reference network for pint with all projects removed
# maybe worth to merge with toot rule, if they turn out to be very similar
rule prepare_pint_reference:
    input:
        network=rules.simplify_sb_network.output.network,
        transmission_projects=rules.clean_projects.output.transmission_projects,
        storage_projects=rules.clean_projects.output.storage_projects,
    output:
        network=resources("cba/pint/networks/reference_{planning_horizons}.nc"),
    script:
        "../../scripts/cba/prepare_pint_reference.py"


# remove the single project {cba_project} from the toot reference network
# currently this can be either a trans123 or a stor123 project
rule prepare_toot_project:
    input:
        network=rules.prepare_toot_reference.output.network,
        transmission_projects=rules.clean_projects.output.transmission_projects,
        storage_projects=rules.clean_projects.output.storage_projects,
    output:
        network=resources(
            "cba/toot/networks/project_{cba_project}_{planning_horizons}.nc"
        ),
    script:
        "../../scripts/cba/prepare_toot_project.py"


# add the single project {cba_project} to the pint reference network
# currently this can be either a trans123 or a stor123 project
rule prepare_pint_project:
    input:
        network=rules.prepare_pint_reference.output.network,
        transmission_projects=rules.clean_projects.output.transmission_projects,
        storage_projects=rules.clean_projects.output.storage_projects,
    output:
        network=resources(
            "cba/pint/networks/project_{cba_project}_{planning_horizons}.nc"
        ),
    script:
        "../../scripts/cba/prepare_pint_project.py"


# solve any of the prepared networks, ie a reference or a project network
# should reuse/import functions from solve_network.py
rule solve_cba_network:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        time_resolution=config_provider("clustering", "temporal", "resolution_sector"),
    input:
        network=resources("cba/{cba_method}/networks/{name}_{planning_horizons}.nc"),
    output:
        network=resources("cba/{cba_method}/postnetworks/{name}_{planning_horizons}.nc"),
    script:
        "../../scripts/cba/solve_cba_network.py"


# compute all metrics for a single pint or toot project comparing reference and project solution
rule make_indicators:
    input:
        reference=resources(
            "cba/{cba_method}/postnetworks/reference_{planning_horizons}.nc"
        ),
        project=resources(
            "cba/{cba_method}/postnetworks/project_{cba_project}_{planning_horizons}.nc"
        ),
    output:
        indicators=RESULTS
        + "cba/{cba_method}/project_{cba_project}_{planning_horizons}.csv",
    script:
        "../../scripts/cba/make_indicators.py"


def input_indicators(w):
    """
    List all indicators csv
    """
    transmission_projects = pd.read_csv(
        checkpoints.clean_projects.get(**w).output.transmission_projects
    )
    storage_projects = pd.read_csv(
        checkpoints.clean_projects.get(**w).output.storage_projects
    )
    planning_horizons = config_provider("scenario", "planning_horizons")(w)

    cba_projects = [
        f"t{pid}" for pid in transmission_projects["project_id"].unique()
    ] + [f"s{pid}" for pid in storage_projects["project_id"].unique()]

    return expand(
        rules.make_indicators.output.indicators,
        cba_project=cba_projects,
        cba_method=["toot", "pint"],  # maybe promote to config
        planning_horizons=planning_horizons,
        allow_missing=True,
    )


# collect the indicators for all transmission_projects into a single overview csv
rule collect_indicators:
    input:
        indicators=input_indicators,
    output:
        indicators=RESULTS + "cba/indicators.csv",
    script:
        "../../scripts/cba/collect_indicators.py"


# pseudo-rule, to run enable running cba with snakemake cba --configfile config/config.tyndp.yaml
rule cba:
    input:
        lambda w: expand(
            rules.collect_indicators.output.indicators,
            run=config_provider("run", "name")(w),
        ),
