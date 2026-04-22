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
    run="(?!None)[-a-zA-Z0-9]+",  # Disallow None as a run wildcard


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


def _effective_horizon(h, warn_fn=None, msg=None):
    if h not in [2030, 2040]:
        if warn_fn:
            warn_fn(msg or "Using 2040 for unsupported planning horizon %s.", h)
        return 2040
    return h


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
        buses=rules.retrieve_tyndp.output.nodes,
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


def input_sb_network(w, run=None):
    """Give correct SB input network path for CBA.

    If config["cba"]["cba_scenario_input"]["use_presolved"] is True,
    use the path to the presolved SB network based on the planning horizon.

    Otherwise, returns the path to the regular SB input network
    based on the scenario wildcards and foresight."""
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
        horizon = _effective_horizon(
            int(w.planning_horizons),
            warn_fn=logger.warning,
            msg=(
                "Presolved SB networks are only available for 2030 and 2040. "
                "Falling back to 2040 for CBA planning horizon %s."
            ),
        )
        return presolved_sb_network_path(w, horizon)

    expanded_wildcards = {
        "clusters": clusters,
        "opts": opts,
        "sector_opts": sector_opts,
    }
    if run is not None:
        expanded_wildcards["run"] = run

    match config_provider("foresight")(w):
        case "perfect":
            expanded_wildcards["planning_horizons"] = "all"
        case "myopic":
            expanded_wildcards["planning_horizons"] = _effective_horizon(
                int(w.planning_horizons),
                warn_fn=logger.warning,
                msg=(
                    "CBA planning horizon %s is not supported for SB inputs. "
                    "Using 2040 inputs instead."
                ),
            )
            # converts the same value to a string so fill_wildcards() can safely call .replace().
            expanded_wildcards["planning_horizons"] = str(
                expanded_wildcards["planning_horizons"]
            )
        case _:
            raise ValueError('config["foresight"] must be one of "perfect" or "myopic"')

    return fill_wildcards(
        RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
        **expanded_wildcards,
    )


def cba_sb_run(w):
    """Determine run name to use as the SB input network for CBA.

    Logic is as follows:
    - If config["cba"]["sb_scenario"] is set, use that run (allows overriding for testing or special cases)
    - Otherwise, if the current run has a "run" wildcard, use that (allows regular scenarios to use their own run as SB input).
    """
    sb_run = config_provider("cba", "sb_scenario", default=None)(w)
    if sb_run:
        return sb_run
    current_run = w.get("run", config_provider("run", "name")(w))
    if isinstance(current_run, list):
        current_run = current_run[0] if current_run else ""
    if current_run:
        return current_run
    return config_provider("tyndp_scenario")(w)


def cba_is_first_horizon(w):
    """Return True if the current horizon is the first CBA horizon.

    Within the CBA workflow, for a multi-climate year run, for the first
    planning horizon, the network output of ``add_existing_baseyear`` is
    used as the climate-year source network, while for later horizons,
    the network output of ``prepare_sector_network`` is used.

    This function checks if the current horizon is the first one to
    determine which source network to use.
    """
    planning_horizons = config_provider("scenario", "planning_horizons")(w)
    if not planning_horizons:
        raise ValueError("scenario.planning_horizons must contain at least one value.")
    return int(w.planning_horizons) == int(planning_horizons[0])


def cba_climate_source_network(w):
    """
    Return the climate-year source network used to overwrite temporal inputs for CBA.

    First horizon:
    - use the output of ``add_existing_baseyear`` because it already contains the
      target climate year's pre-solve state without depending on a previous
      climate-year solve.

    Later horizons:
    - use the output of ``prepare_sector_network``. A dedicated CBA-side
      carry-over/mapping step is still needed there, so the preparation script
      handles that case explicitly.
    """
    scenario = config_provider("scenario")(w)
    (clusters,) = scenario["clusters"]
    (opts,) = scenario["opts"]
    (sector_opts,) = scenario["sector_opts"]
    run = w.get("run", config_provider("run", "name")(w))
    if isinstance(run, list):
        run = run[0] if run else ""
    if cba_is_first_horizon(w):
        return fill_wildcards(
            resources(
                "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}_brownfield.nc"
            ),
            clusters=clusters,
            opts=opts,
            sector_opts=sector_opts,
            planning_horizons=str(w.planning_horizons),
            run=run,
        )
    return fill_wildcards(
        resources(
            "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc"
        ),
        clusters=clusters,
        opts=opts,
        sector_opts=sector_opts,
        planning_horizons=str(w.planning_horizons),
        run=run,
    )


# Prepare the correct SB input network for CBA, depending on the horizon and whether presolved networks are used
rule prepare_cba_sb_input_network:
    input:
        solved_network=lambda w: input_sb_network(w, run=cba_sb_run(w)),
        climate_source_network=cba_climate_source_network,
    params:
        is_first_horizon=cba_is_first_horizon,
    output:
        network=resources("cba/networks/sb_input_{planning_horizons}.nc"),
    log:
        logs("cba/prepare_cba_sb_input_network_{planning_horizons}.log"),
    script:
        scripts("cba/prepare_cba_sb_input_network.py")


# Simplify scenario building network for CBA
# Fixes capacities, adds hurdle costs, extends primary fuel sources, disables volume limits
rule simplify_sb_network:
    params:
        tyndp_conventional_carriers=config_provider(
            "electricity", "tyndp_conventional_carriers"
        ),
        hurdle_costs=config_provider("cba", "hurdle_costs"),
    input:
        network=rules.prepare_cba_sb_input_network.output.network,
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
        buses=rules.build_tyndp_network.output.substations_geojson,
    output:
        corrections=resources("cba/reference_sb_to_cba_{planning_horizons}.csv"),
    log:
        logs("cba/fix_reference_sb_to_cba_{planning_horizons}.log"),
    benchmark:
        benchmarks("performances/cba/fix_reference_sb_to_cba_{planning_horizons}")
    script:
        scripts("cba/fix_reference_sb_to_cba.py")


# Build reference network with all TOOT projects included
# Ensures MSV extraction and rolling horizon use the same baseline
rule prepare_reference:
    params:
        hurdle_costs=config_provider("cba", "hurdle_costs"),
        patch_sb_with_annexe=config_provider(
            "tyndp_investment_candidates", "patch_sb_with_annexe"
        ),
    input:
        network=rules.simplify_sb_network.output.network,
        transmission_projects=rules.clean_projects.output.transmission_projects,
        storage_projects=rules.clean_projects.output.storage_projects,
        corrections=rules.fix_reference_sb_to_cba.output.corrections,
        costs=resources("costs_{planning_horizons}_processed.csv"),
    output:
        network=resources("cba/networks/reference_{planning_horizons}.nc"),
    script:
        scripts("cba/prepare_reference.py")


# Generate snapshot weightings for MSV extraction temporal aggregation
rule build_msv_snapshot_weightings:
    params:
        msv_resolution=config_provider("cba", "msv_extraction", "resolution"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
    input:
        network=rules.prepare_reference.output.network,
    output:
        snapshot_weightings=resources(
            "cba/msv_snapshot_weightings_{planning_horizons}.csv"
        ),
    script:
        "../scripts/cba/build_msv_snapshot_weightings.py"


def input_msv_snapshot_weightings(w):
    """Return snapshot weightings file only if MSV resolution requires it."""
    resolution = config_provider("cba", "msv_extraction", "resolution")(w)
    if resolution and "h" in str(resolution).lower():
        return rules.build_msv_snapshot_weightings.output.snapshot_weightings
    return []


# Extract marginal storage values via perfect foresight solve (full year with cyclicity enabled)
rule solve_cba_msv_extraction:
    params:
        solving=config_provider("solving"),
        cba_solving=config_provider("cba", "msv_extraction", "solving"),
        msv_resolution=config_provider("cba", "msv_extraction", "resolution"),
        cyclic_carriers=config_provider("cba", "storage", "cyclic_carriers"),
    input:
        network=rules.prepare_reference.output.network,
        snapshot_weightings=input_msv_snapshot_weightings,
    output:
        network=resources("cba/networks/msv_{planning_horizons}.nc"),
    log:
        solver=RESULTS + "logs/cba/msv/{planning_horizons}_solver.log",
        memory=RESULTS + "logs/cba/msv/{planning_horizons}_memory.log",
        python=RESULTS + "logs/cba/msv/{planning_horizons}_python.log",
    threads: 1
    script:
        "../scripts/cba/solve_cba_msv_extraction.py"


# Prepare network for rolling horizon: disable seasonal cyclicity, apply marginal storage value
rule prepare_rolling_horizon:
    params:
        cyclic_carriers=config_provider("cba", "storage", "cyclic_carriers"),
        soc_boundary_carriers=config_provider("cba", "storage", "soc_boundary_carriers"),
        msv_resample_method=config_provider("cba", "msv_extraction", "resample_method"),
    input:
        network=rules.prepare_reference.output.network,
        network_msv=rules.solve_cba_msv_extraction.output.network,
    output:
        network=resources("cba/networks/rl_{planning_horizons}.nc"),
    script:
        scripts("cba/prepare_rolling_horizon.py")


# add or remove the cba project based on assigned method
rule prepare_project:
    params:
        hurdle_costs=config_provider("cba", "hurdle_costs"),
        cyclic_carriers=config_provider("cba", "storage", "cyclic_carriers"),
        soc_boundary_carriers=config_provider("cba", "storage", "soc_boundary_carriers"),
    input:
        network=rules.prepare_rolling_horizon.output.network,
        network_msv=rules.solve_cba_msv_extraction.output.network,
        transmission_projects=rules.clean_projects.output.transmission_projects,
        storage_projects=rules.clean_projects.output.storage_projects,
        methods=rules.clean_projects.output.methods,
        costs=resources("costs_{planning_horizons}_processed.csv"),
    output:
        network=resources("cba/networks/project_{cba_project}_{planning_horizons}.nc"),
    script:
        scripts("cba/prepare_project.py")


# Solve reference network with rolling horizon (MSV already applied)
rule solve_cba_reference_network:
    params:
        solving=config_provider("solving"),
        cba_solving=config_provider("cba", "solving"),
        foresight=config_provider("foresight"),
        time_resolution=config_provider("clustering", "temporal", "resolution_sector"),
        custom_extra_functionality=None,
    input:
        network=rules.prepare_rolling_horizon.output.network,
    output:
        network=RESULTS + "cba/networks/reference_{planning_horizons}.nc",
    log:
        solver=RESULTS + "logs/cba/reference/reference_{planning_horizons}_solver.log",
        memory=RESULTS + "logs/cba/reference/reference_{planning_horizons}_memory.log",
        python=RESULTS + "logs/cba/reference/reference_{planning_horizons}_python.log",
    threads: 1
    script:
        scripts("cba/solve_cba_network.py")


# Solve TOOT/PINT project network with rolling horizon
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
        solver=RESULTS
        + "logs/cba/projects/project_{cba_project}_{planning_horizons}_solver.log",
        memory=RESULTS
        + "logs/cba/projects/project_{cba_project}_{planning_horizons}_memory.log",
        python=RESULTS
        + "logs/cba/projects/project_{cba_project}_{planning_horizons}_python.log",
    threads: 1
    script:
        scripts("cba/solve_cba_network.py")


# Compute CBA indicators comparing reference and project networks
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

    Works for collection scenarios and regular scenarios.
    """
    run = w.get("run", config_provider("run", "name")(w))
    projects = pd.read_csv(checkpoints.clean_projects.get(run=run).output.methods)
    horizon = _effective_horizon(
        int(w.planning_horizons),
        warn_fn=logger.warning,
        msg=(
            "CBA methods are only available for 2030 or 2040. "
            "Using 2040 for planning horizon %s."
        ),
    )
    if "planning_horizon" in projects.columns:
        projects = projects.loc[projects["planning_horizon"] == horizon]
    cba_projects = [f"t{pid}" for pid in projects["project_id"].unique()]

    # Collection scenarios look for results within regular scenarios
    runs = config_provider("cba", "scenarios", default="{run}")
    # NOTE: project_specs filtering happens on the collection scenario (and does not descend)
    project_specs = config_provider("cba", "projects")(w)

    return expand(
        rules.make_indicators.output.indicators,
        cba_project=filter_projects_by_specs(cba_projects, project_specs),
        run=runs,
        allow_missing=True,
    )


# Collect indicators for all projects into overview CSV
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


# rule plot_all_cba_benchmark:
#     input:
#         indicators=rules.collect_indicators.output.indicators,
#     output:
#         plot_dir=directory(RESULTS + "cba/validation_{planning_horizons}"),
#     script:
#         scripts("cba/plot_benchmark_indicators.py")


rule plot_weather_benchmark:
    input:
        # indicators=rules.collect_indicators.output.indicators,
        indicators=rules.make_indicators.output.indicators,
    output:
        plot_file=RESULTS
        + "cba/ensemble_plots/ensemble_{cba_project}_{planning_horizons}.png",
    script:
        "../scripts/cba/plot_benchmark_indicators.py"


rule average_indicators_per_project_and_planning_horizon:
    input:
        indicators=lambda w: expand(
            rules.make_indicators.output.indicators,
            run=config_provider("cba", "scenarios")(w),
            allow_missing=True,
        ),
    output:
        indicators=RESULTS
        + "cba/ensemble_indicators/ensemble_indicators_{cba_project}_{planning_horizons}.csv",
    script:
        "../scripts/cba/average_indicators.py"


rule summarize_indicators_per_project:
    input:
        indicators=lambda w: expand(
            rules.average_indicators_per_project_and_planning_horizon.output.indicators,
            transmission_projects=rules.clean_projects.output.transmission_projects,
            planning_horizons=config["cba"]["planning_horizons"],
            allow_missing=True,
        ),
    output:
        plot_file=RESULTS + "cba/ensemble_plots/ensemble_{cba_project}_all_horizons.png",
    script:
        "../scripts/cba/summarize_indicators.py"


rule summarize_all_indicators:
    input:
        indicators=lambda w: expand(
            rules.plot_weather_benchmark.input.indicators,
            transmission_projects=rules.clean_projects.output.transmission_projects,
            planning_horizons=config["cba"]["planning_horizons"],
            cba_project=cba_projects(w),
            run=config_provider("cba", "scenarios")(w),
        ),
    output:
        plot_file=RESULTS + "cba/ensemble_plots/ensemble_all.png",
    script:
        "../scripts/cba/summarize_all.py"


# pseudo-rule, to run enable running cba with snakemake cba --configfile config/config.tyndp.yaml
def cba_collection_scenarios(w):
    """
    Return cba collection scenarios, ie. the meta scenarios with cba: scenarios config
    """
    names = config["run"]["name"]
    if isinstance(names, str):
        names = [names]
    scenarios = []
    # fall back to the raw run.name if it isn’t found in the scenarios file
    for name in names:
        try:
            scn = scenario_config(name)
        except KeyError:
            scenarios.append(name)
            continue
        if scn.get("cba", {}).get("scenarios") is not None:
            scenarios.append(name)
    return scenarios


def cba_scenarios(w):
    """
    Return cba: scenarios of a cba collection scenario
    """
    run = w.get("run", config_provider("run", "name")(w))
    if isinstance(run, list):
        run = run[0] if run else ""
    try:
        scn = scenario_config(run)
    except KeyError:
        return [run] if run else []
    return scn.get("cba", {}).get("scenarios", [run] if run else [])


def cba_projects(w):
    """
    List all indicators csv
    """
    # run = config_provider("run", "name")(w),
    run = w.get("run", config_provider("run", "name")(w))
    if isinstance(run, list):
        run = run[0] if run else ""
    if "cy" not in run and run:
        run = f"{run}-cy2009"

    projects = pd.read_csv(checkpoints.clean_projects.get(run=run).output.methods)
    cba_projects = [f"t{pid}" for pid in projects["project_id"].unique()]
    project_specs = config_provider("cba", "projects")(w)
    cba_project = filter_projects_by_specs(cba_projects, project_specs)

    return expand(cba_project)


# collect files to be stored in the scenario directory, e.g., NT-cy1995
rule collect_cba_scenario:
    input:
        lambda w: expand(
            rules.plot_weather_benchmark.output.plot_file,
            planning_horizons=config_provider("cba", "planning_horizons")(w),
            cba_project=cba_projects(w),
            run=cba_scenarios(w),
        ),
        lambda w: expand(
            rules.plot_indicators.output.plot_dir,
            planning_horizons=config_provider("cba", "planning_horizons")(w),
            run=cba_scenarios(w),
        ),
        lambda w: expand(
            rules.plot_cba_benchmark.output.plot_file,
            planning_horizons=config_provider("cba", "planning_horizons")(w),
            cba_project=cba_projects(w),
            run=cba_scenarios(w),
        ),
    output:
        touch(RESULTS + "cba/all_scenarios.txt"),


# collect files to be stored in the scenario collection directory, e.g., NT-cyears
rule cba:
    input:
        lambda w: expand(
            rules.average_indicators_per_project_and_planning_horizon.output.indicators,
            planning_horizons=config["cba"]["planning_horizons"],
            cba_project=cba_projects(w),
            run=cba_collection_scenarios(w),
        ),
        lambda w: expand(
            rules.summarize_indicators_per_project.output.plot_file,
            cba_project=cba_projects(w),
            run=cba_collection_scenarios(w),
        ),
        lambda w: expand(
            rules.summarize_all_indicators.output.plot_file,
            planning_horizons=config["cba"]["planning_horizons"],
            cba_project=cba_projects(w),
            run=cba_collection_scenarios(w),
        ),
        # lambda w: expand(
        #     rules.plot_all_cba_benchmark.output.plot_dir,
        #     planning_horizons=config["cba"]["planning_horizons"],
        #     run=cba_collection_scenarios(w),
        # ),
        # collect files to be stored in the scenario directory, e.g., NT-cy1995
        lambda w: expand(
            rules.collect_cba_scenario.output[0],
            run=cba_collection_scenarios(w),
        ),


# collect rules
rule prepare_references:
    input:
        lambda w: expand(
            resources("cba/networks/reference_{planning_horizons}.nc"),
            **config["scenario"],
            run=config["run"]["name"],
        ),
