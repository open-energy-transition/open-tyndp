# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
#

import logging
import os
import shutil
from pathlib import Path
from zipfile import ZipFile

import pandas as pd

from scripts.cba._helpers import filter_projects_by_specs
from scripts._helpers import fill_wildcards
from shutil import unpack_archive, copy2

logger = logging.getLogger(__name__)


wildcard_constraints:
    cba_project=r"(s|t)\d+",


if (CBA_PROJECTS_DATASET := dataset_version("tyndp_cba_projects"))["source"] in [
    "archive"
]:

    rule retrieve_tyndp_cba_projects:
        params:
            source="CBA project explorer",
        input:
            zip_file=storage(CBA_PROJECTS_DATASET["url"]),
        output:
            dir=directory(CBA_PROJECTS_DATASET["folder"]),
        log:
            "logs/retrieve_tyndp_cba_projects.log",
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
            file=f"{CBA_NON_CO2_DATASET["folder"]}/a.3_non-co2-emissions.csv",
        log:
            "logs/retrieve_tyndp_cba_non_co2_emissions.log",
        run:
            copy2(input["file"], output["file"])


if (CBA_GUIDELINES_DATASET := dataset_version("cba_guidelines_reference_projects"))[
    "source"
] in ["archive"]:

    rule retreive_cba_guidelines_reference_projects:
        input:
            file=storage(CBA_GUIDELINES_DATASET["url"]),
        output:
            file=f"{CBA_GUIDELINES_DATASET['folder']}/table_B1_CBA_Implementations_Guidelines_TYNDP2024.csv",
        log:
            "logs/retreive_cba_guidelines_reference_projects.log",
        run:
            copy2(input["file"], output["file"])


def presolved_sb_network_path(w, horizon=None):
    sb_version = config_provider(
        "cba", "cba_scenario_input", "sb_version", default="latest"
    )(w)
    target_horizon = horizon if horizon is not None else w.planning_horizons
    return RESULTS + f"networks/presolved-{sb_version}/base_s_all___{target_horizon}.nc"


if config.get("cba", {}).get("cba_scenario_input", {}).get("use_presolved", False):
    if (SB_SOLVED_NETWORKS_DATASET := dataset_version("open_tyndp_prelim"))[
        "source"
    ] in ["archive"]:

        rule retrieve_presolved_sb_networks:
            input:
                zip_file=storage(SB_SOLVED_NETWORKS_DATASET["url"]),
            output:
                network=RESULTS
                + f"networks/presolved-{config['cba']['cba_scenario_input']['sb_version']}/base_s_all___{{planning_horizons}}.nc",
            log:
                logs("retrieve_presolved_sb_networks_{planning_horizons}.log"),
            run:
                target_suffix = (
                    f"networks/base_s_all___{wildcards.planning_horizons}.nc"
                )
                with ZipFile(input["zip_file"], "r") as zf:
                    matches = [
                        m for m in zf.namelist() if m.endswith(target_suffix)
                    ]
                    if not matches:
                        raise ValueError(
                            f"Could not find '{target_suffix}' in {input['zip_file']}."
                        )
                    member = matches[0]
                    out_path = Path(output["network"])
                    out_path.parent.mkdir(parents=True, exist_ok=True)
                    with zf.open(member) as src, out_path.open("wb") as dst:
                        shutil.copyfileobj(src, dst)



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
        guidelines=rules.retreive_cba_guidelines_reference_projects.output.file,
    output:
        # TODO: The toot_projects and pint_projects outputs are likely only
        # transmission projects (no storage). In order to confirm, we should check 
        # if Table B.1 from the guidelines (table_B1_CBA_Implementations_Guidelines_TYNDP2024.csv) 
        # contains only transmission or also storage projects.
        transmission_projects=resources("cba/transmission_projects.csv"),
        storage_projects=resources("cba/storage_projects.csv"),
        methods=resources("cba/cba_project_methods.csv"),
    script:
        scripts("cba/clean_projects.py")


rule clean_tyndp_indicators:
    input:
        dir=rules.retrieve_tyndp_cba_projects.output.dir,
    output:
        indicators=resources("cba/tyndp_indicators.csv"),
        readme=resources("cba/tyndp_indicators_name_unit.csv"),
    script:
        scripts("cba/clean_tyndp_indicators.py")


def input_sb_network(w):
    scenario = config_provider("scenario")(w)
    (clusters,) = scenario["clusters"]
    (opts,) = scenario["opts"]
    (sector_opts,) = scenario["sector_opts"]

    if config_provider("cba", "cba_scenario_input", "use_presolved", default=False)(w):
        scenario_name = config_provider("tyndp_scenario")(w)
        if scenario_name != "NT":
            raise ValueError(
                "Presolved SB networks are only currently available for the NT scenario."
            )
        sb_version = config_provider(
            "cba", "cba_scenario_input", "sb_version", default="latest"
        )(w)
        # If sb_version is not "latest", raise error (for now) until we add functionality to handle gathering different Zenodo versions
        if sb_version != "latest":
            raise ValueError(
                "Only cba.cba_scenario_input.sb_version='latest' is supported for presolved SB runs at the moment."
            )
        # Check that options match the presolved network naming convention
        if clusters != "all" or opts != "" or sector_opts != "":
            raise ValueError(
                "Presolved SB runs require scenario.clusters=['all'], "
                "scenario.opts=[''], and scenario.sector_opts=[''] to match "
                "the Zenodo network naming (base_s_all___{planning_horizons}.nc)."
            )
        horizon = int(w.planning_horizons)
        # If CBA planning horizon not in [2030, 2040] (such as horizon == 2035), use the 2040 SB network
        if horizon not in [2030, 2040]:
            logger.warning(
                "Presolved SB networks are only available for 2030 and 2040. "
                "Falling back to 2040 for CBA planning horizon %s.",
                w.planning_horizons,
            )
            horizon = 2040
        return presolved_sb_network_path(w, horizon)

    expanded_wildcards = {
        "clusters": clusters,
        "opts": opts,
        "sector_opts": sector_opts,
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
        scripts("cba/simplify_sb_network.py")


# build reference corrections between SB investments and CBA guidelines


def get_elec_project_build_years(w):
    return config_provider("tyndp_investment_candidates", "elec_projects")(w).get(
        int(w.planning_horizons)
    )


rule fix_reference_sb_to_cba:
    params:
        build_years=get_elec_project_build_years,
    input:
        invest_grid=rules.retrieve_tyndp.output.invest_grid,
        guidelines=rules.retreive_cba_guidelines_reference_projects.output.file,
        transmission_projects=rules.clean_projects.output.transmission_projects,
    output:
        corrections_csv=resources("cba/reference_sb_to_cba_{planning_horizons}.csv"),
    log:
        logs("cba/fix_reference_sb_to_cba_{planning_horizons}.log"),
    benchmark:
        benchmarks("cba/fix_reference_sb_to_cba_{planning_horizons}")
    script:
        scripts("cba/fix_reference_sb_to_cba.py")


# build the reference network
rule prepare_reference:
    params:
        hurdle_costs=config_provider("cba", "hurdle_costs"),
    input:
        network=rules.simplify_sb_network.output.network,
        transmission_projects=rules.clean_projects.output.transmission_projects,
        storage_projects=rules.clean_projects.output.storage_projects,
        corrections=rules.fix_reference_sb_to_cba.output.corrections_csv,
    output:
        network=resources("cba/networks/reference_{planning_horizons}.nc"),
    script:
        scripts("cba/prepare_reference.py")


# add or remove the cba project based on assigned method
rule prepare_project:
    params:
        hurdle_costs=config_provider("cba", "hurdle_costs"),
    input:
        network=rules.prepare_reference.output.network,
        transmission_projects=rules.clean_projects.output.transmission_projects,
        storage_projects=rules.clean_projects.output.storage_projects,
        methods=rules.clean_projects.output.methods,
    output:
        network=resources("cba/networks/project_{cba_project}_{planning_horizons}.nc"),
    script:
        scripts("cba/prepare_project.py")


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
        scripts("cba/solve_cba_network.py")


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
        network=resources("cba/networks/project_{cba_project}_{planning_horizons}.nc"),
    output:
        network=RESULTS + "cba/networks/project_{cba_project}_{planning_horizons}.nc",
    log:
        solver=logs(
            "cba/solve_cba_network/project_{cba_project}_{planning_horizons}_solver.log"
        ),
        memory=logs(
            "cba/solve_cba_network/project_{cba_project}_{planning_horizons}_memory.log"
        ),
        python=logs(
            "cba/solve_cba_network/project_{cba_project}_{planning_horizons}_python.log"
        ),
    threads: 1
    script:
        scripts("cba/solve_cba_network.py")


# compute all metrics for a single pint or toot project comparing reference and project solution
rule make_indicators:
    input:
        reference=RESULTS + "cba/networks/reference_{planning_horizons}.nc",
        project=RESULTS + "cba/networks/project_{cba_project}_{planning_horizons}.nc",
        non_co2_emissions=rules.retrieve_tyndp_cba_non_co2_emissions.output.file,
        benchmark=rules.clean_tyndp_indicators.output.indicators,
        methods=rules.clean_projects.output.methods,
    output:
        indicators=RESULTS + "cba/project_{cba_project}_{planning_horizons}.csv",
    script:
        scripts("cba/make_indicators.py")


def input_indicators(w):
    """
    List all indicators csv
    """
    run = w.get("run", config_provider("run", "name")(w))
    projects = pd.read_csv(checkpoints.clean_projects.get(run=run).output.methods)
    if "planning_horizon" in projects.columns:
        projects = projects.loc[
            projects["planning_horizon"] == int(w.planning_horizons)
        ]
    cba_projects = [f"t{pid}" for pid in projects["project_id"].unique()]

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
        indicators=RESULTS + "cba/indicators_{planning_horizons}.csv",
    script:
        scripts("cba/collect_indicators.py")


rule plot_indicators:
    params:
        plotting=config_provider("plotting"),
    input:
        indicators=rules.collect_indicators.output.indicators,
        transmission_projects=rules.clean_projects.output.transmission_projects,
    output:
        plot_dir=directory(RESULTS + "cba/plots_{planning_horizons}"),
    script:
        scripts("cba/plot_indicators.py")


rule plot_cba_benchmark:
    input:
        indicators=RESULTS + "cba/project_{cba_project}_{planning_horizons}.csv",
    output:
        plot_file=RESULTS
        + "cba/validation_{planning_horizons}/project_{cba_project}_{planning_horizons}.png",
    script:
        scripts("cba/plot_benchmark_indicators.py")


rule plot_all_cba_benchmark:
    input:
        indicators=rules.collect_indicators.output.indicators,
    output:
        plot_dir=directory(RESULTS + "cba/validation_{planning_horizons}"),
    script:
        scripts("cba/plot_benchmark_indicators.py")


# pseudo-rule, to run enable running cba with snakemake cba --configfile config/config.tyndp.yaml
rule cba:
    input:
        lambda w: expand(
            rules.collect_indicators.output.indicators,
            planning_horizons=config_provider("cba", "planning_horizons")(w),
            run=config_provider("run", "name")(w),
        ),
        lambda w: expand(
            rules.plot_indicators.output.plot_dir,
            planning_horizons=config_provider("cba", "planning_horizons")(w),
            run=config_provider("run", "name")(w),
        ),
        lambda w: expand(
            rules.plot_all_cba_benchmark.output.plot_dir,
            planning_horizons=config_provider("cba", "planning_horizons")(w),
            run=config_provider("run", "name")(w),
        ),


# collect rules
rule prepare_references:
    input:
        lambda w: expand(
            resources("cba/networks/reference_{planning_horizons}.nc"),
            **config["scenario"],
            run=config["run"]["name"],
        ),
