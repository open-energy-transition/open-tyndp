# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import contextlib
import copy
import hashlib
import logging
import os
import re
import time
from bisect import bisect_right
from collections.abc import Callable
from functools import partial, wraps
from pathlib import Path
from tempfile import NamedTemporaryFile

import atlite
import fiona
import git
import numpy as np
import pandas as pd
import pypsa
import pytz
import requests
import xarray as xr
import yaml
from snakemake.utils import update_config
from tenacity import (
    retry as tenacity_retry,
)
from tenacity import (
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)
from tqdm import tqdm

logger = logging.getLogger(__name__)

REGION_COLS = ["geometry", "name", "x", "y", "country"]

SCENARIO_DICT = {
    "Distributed Energy": "DE",
    "Global Ambition": "GA",
    r"National Trends\s*\+": "NT",
    r"NT\s*\+": "NT",
    "National Trends": "NT",
}

ENERGY_UNITS = {"TWh", "GWh", "MWh", "kWh"}
POWER_UNITS = {"GW", "MW", "kW"}

PYPSA_V1 = bool(re.match(r"^1\.\d", pypsa.__version__))


def get_scenarios(run):
    scenario_config = run.get("scenarios", {})
    if run["name"] and scenario_config.get("enable"):
        fn = Path(scenario_config["file"])
        if fn.exists():
            scenarios = yaml.safe_load(fn.read_text())
            if run["name"] == "all":
                run["name"] = list(scenarios.keys())
            return scenarios
    return {}


def get_rdir(run):
    scenario_config = run.get("scenarios", {})
    if run["name"] and scenario_config.get("enable"):
        RDIR = "{run}/"
    elif run["name"]:
        RDIR = run["name"] + "/"
    else:
        RDIR = ""

    prefix = run.get("prefix", "")
    if prefix:
        RDIR = f"{prefix}/{RDIR}"

    return RDIR


def get_run_path(fn, dir, rdir, shared_resources, exclude_from_shared):
    """
    Dynamically provide paths based on shared resources and filename.

    Use this function for snakemake rule inputs or outputs that should be
    optionally shared across runs or created individually for each run.

    Parameters
    ----------
    fn : str
        The filename for the path to be generated.
    dir : str
        The base directory.
    rdir : str
        Relative directory for non-shared resources.
    shared_resources : str or bool
        Specifies which resources should be shared.
        - If string is "base", special handling for shared "base" resources (see notes).
        - If random string other than "base", this folder is used instead of the `rdir` keyword.
        - If boolean, directly specifies if the resource is shared.
    exclude_from_shared: list
        List of filenames to exclude from shared resources. Only relevant if shared_resources is "base".

    Returns
    -------
    str
        Full path where the resource should be stored.

    Notes
    -----
    Special case for "base" allows no wildcards other than "technology", "year"
    and "scope" and excludes filenames starting with "networks/elec" or
    "add_electricity". All other resources are shared.
    """
    if shared_resources == "base":
        pattern = r"\{([^{}]+)\}"
        existing_wildcards = set(re.findall(pattern, fn))
        irrelevant_wildcards = {"technology", "year", "scope", "kind"}
        no_relevant_wildcards = not existing_wildcards - irrelevant_wildcards
        not_shared_rule = (
            not fn.endswith("elec.nc")
            and not fn.startswith("add_electricity")
            and not any(fn.startswith(ex) for ex in exclude_from_shared)
        )
        is_shared = no_relevant_wildcards and not_shared_rule
        shared_files = (
            "networks/base_s_{clusters}.nc",
            "regions_onshore_base_s_{clusters}.geojson",
            "regions_offshore_base_s_{clusters}.geojson",
            "busmap_base_s_{clusters}.csv",
            "linemap_base_s_{clusters}.csv",
            "cluster_network_base_s_{clusters}",
            "profile_{clusters}_",
            "build_renewable_profile_{clusters}",
            "regions_by_class_{clusters}",
            "availability_matrix_",
            "determine_availability_matrix_",
            "solar_thermal",
        )
        if any(prefix in fn for prefix in shared_files) or is_shared:
            is_shared = True
        rdir = "" if is_shared else rdir
    elif isinstance(shared_resources, str):
        rdir = shared_resources + "/"
    elif isinstance(shared_resources, bool):
        rdir = "" if shared_resources else rdir
    else:
        raise ValueError(
            "shared_resources must be a boolean, str, or 'base' for special handling."
        )

    return f"{dir}{rdir}{fn}"


def path_provider(dir, rdir, shared_resources, exclude_from_shared):
    """
    Returns a partial function that dynamically provides paths based on shared
    resources and the filename.

    Returns
    -------
    partial function
        A partial function that takes a filename as input and
        returns the path to the file based on the shared_resources parameter.
    """
    return partial(
        get_run_path,
        dir=dir,
        rdir=rdir,
        shared_resources=shared_resources,
        exclude_from_shared=exclude_from_shared,
    )


def get_shadow(run):
    """
    Returns 'shallow' or None depending on the user setting.
    """
    shadow_config = run.get("use_shadow_directory", True)
    if shadow_config:
        return "shallow"
    return None


def get_opt(opts, expr, flags=None):
    """
    Return the first option matching the regular expression.

    The regular expression is case-insensitive by default.
    """
    if flags is None:
        flags = re.IGNORECASE
    for o in opts:
        match = re.match(expr, o, flags=flags)
        if match:
            return match.group(0)
    return None


def find_opt(opts, expr):
    """
    Return if available the float after the expression.
    """
    for o in opts:
        if expr in o:
            m = re.findall(r"m?\d+(?:[\.p]\d+)?", o)
            if len(m) > 0:
                return True, float(m[-1].replace("p", ".").replace("m", "-"))
            else:
                return True, None
    return False, None


def fill_wildcards(s: str, **wildcards: str) -> str:
    """
    Fill given (subset of) wildcards into a path with wildcards
    """
    for k, v in wildcards.items():
        if isinstance(v, (list, tuple)):
            assert len(v) == 1, f"Need a single entry for {k}, but found: {v}"
            v = v[0]

        s = s.replace("{" + k + "}", v)
    return s


# Define a context manager to temporarily mute print statements
@contextlib.contextmanager
def mute_print():
    with open(os.devnull, "w") as devnull:
        with contextlib.redirect_stdout(devnull):
            yield


def set_scenario_config(snakemake):
    scenario = snakemake.config["run"].get("scenarios", {})
    if scenario.get("enable") and "run" in snakemake.wildcards.keys():
        try:
            with open(scenario["file"]) as f:
                scenario_config = yaml.safe_load(f)
        except FileNotFoundError:
            # fallback for mock_snakemake
            script_dir = Path(__file__).parent.resolve()
            root_dir = script_dir.parent
            with open(root_dir / scenario["file"]) as f:
                scenario_config = yaml.safe_load(f)
        update_config(snakemake.config, scenario_config[snakemake.wildcards.run])


def configure_logging(snakemake, skip_handlers=False):
    """
    Configure the basic behaviour for the logging module.

    Note: Must only be called once from the __main__ section of a script.

    The setup includes printing log messages to STDERR and to a log file defined
    by either (in priority order): snakemake.log.python, snakemake.log[0] or "logs/{rulename}.log".
    Additional keywords from logging.basicConfig are accepted via the snakemake configuration
    file under snakemake.config.logging.

    Parameters
    ----------
    snakemake : snakemake object
        Your snakemake object containing a snakemake.config and snakemake.log.
    skip_handlers : True | False (default)
        Do (not) skip the default handlers created for redirecting output to STDERR and file.
    """
    import logging
    import sys

    kwargs = snakemake.config.get("logging", dict()).copy()
    kwargs.setdefault("level", "INFO")

    if skip_handlers is False:
        fallback_path = Path(__file__).parent.joinpath(
            "..", "logs", f"{snakemake.rule}.log"
        )
        logfile = snakemake.log.get(
            "python", snakemake.log[0] if snakemake.log else fallback_path
        )
        kwargs.update(
            {
                "handlers": [
                    # Prefer the 'python' log, otherwise take the first log for each
                    # Snakemake rule
                    logging.FileHandler(logfile),
                    logging.StreamHandler(),
                ]
            }
        )
    logging.basicConfig(**kwargs)

    # Setup a function to handle uncaught exceptions and include them with their stacktrace into logfiles
    def handle_exception(exc_type, exc_value, exc_traceback):
        # Log the exception
        logger = logging.getLogger()
        logger.error(
            "Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback)
        )

    sys.excepthook = handle_exception


def update_p_nom_max(n):
    # if extendable carriers (solar/onwind/...) have capacity >= 0,
    # e.g. existing assets from GEM are included to the network,
    # the installed capacity might exceed the expansion limit.
    # Hence, we update the assumptions.

    n.generators.p_nom_max = n.generators[["p_nom_min", "p_nom_max"]].max(1)


def aggregate_p_nom(n):
    return pd.concat(
        [
            n.generators.groupby("carrier").p_nom_opt.sum(),
            n.storage_units.groupby("carrier").p_nom_opt.sum(),
            n.links.groupby("carrier").p_nom_opt.sum(),
            n.loads_t.p.groupby(n.loads.carrier, axis=1).sum().mean(),
        ]
    )


def aggregate_p(n):
    return pd.concat(
        [
            n.generators_t.p.sum().groupby(n.generators.carrier).sum(),
            n.storage_units_t.p.sum().groupby(n.storage_units.carrier).sum(),
            n.stores_t.p.sum().groupby(n.stores.carrier).sum(),
            -n.loads_t.p.sum().groupby(n.loads.carrier).sum(),
        ]
    )


def get(item, investment_year=None):
    """
    Check whether item depends on investment year.
    """
    if not isinstance(item, dict):
        return item
    elif investment_year in item.keys():
        return item[investment_year]
    else:
        logger.warning(
            f"Investment key {investment_year} not found in dictionary {item}."
        )
        keys = sorted(item.keys())
        if investment_year < keys[0]:
            logger.warning(f"Lower than minimum key. Taking minimum key {keys[0]}")
            return item[keys[0]]
        elif investment_year > keys[-1]:
            logger.warning(f"Higher than maximum key. Taking maximum key {keys[0]}")
            return item[keys[-1]]
        else:
            logger.warning(
                "Interpolate linearly between the next lower and next higher year."
            )
            lower_key = max(k for k in keys if k < investment_year)
            higher_key = min(k for k in keys if k > investment_year)
            lower = item[lower_key]
            higher = item[higher_key]
            return lower + (higher - lower) * (investment_year - lower_key) / (
                higher_key - lower_key
            )


def aggregate_e_nom(n):
    return pd.concat(
        [
            (n.storage_units["p_nom_opt"] * n.storage_units["max_hours"])
            .groupby(n.storage_units["carrier"])
            .sum(),
            n.stores["e_nom_opt"].groupby(n.stores.carrier).sum(),
        ]
    )


def aggregate_p_curtailed(n):
    return pd.concat(
        [
            (
                (
                    n.generators_t.p_max_pu.sum().multiply(n.generators.p_nom_opt)
                    - n.generators_t.p.sum()
                )
                .groupby(n.generators.carrier)
                .sum()
            ),
            (
                (n.storage_units_t.inflow.sum() - n.storage_units_t.p.sum())
                .groupby(n.storage_units.carrier)
                .sum()
            ),
        ]
    )


def aggregate_costs(n, flatten=False, opts=None, existing_only=False):
    components = dict(
        Link=("p_nom", "p0"),
        Generator=("p_nom", "p"),
        StorageUnit=("p_nom", "p"),
        Store=("e_nom", "p"),
        Line=("s_nom", None),
        Transformer=("s_nom", None),
    )

    costs = {}
    for c, (p_nom, p_attr) in zip(
        n.iterate_components(components.keys(), skip_empty=False), components.values()
    ):
        if c.df.empty:
            continue
        if not existing_only:
            p_nom += "_opt"
        costs[(c.list_name, "capital")] = (
            (c.df[p_nom] * c.df.capital_cost).groupby(c.df.carrier).sum()
        )
        if p_attr is not None:
            p = c.pnl[p_attr].sum()
            if c.name == "StorageUnit":
                p = p.loc[p > 0]
            costs[(c.list_name, "marginal")] = (
                (p * c.df.marginal_cost).groupby(c.df.carrier).sum()
            )
    costs = pd.concat(costs)

    if flatten:
        assert opts is not None
        conv_techs = opts["conv_techs"]

        costs = costs.reset_index(level=0, drop=True)
        costs = costs["capital"].add(
            costs["marginal"].rename({t: t + " marginal" for t in conv_techs}),
            fill_value=0.0,
        )

    return costs


@tenacity_retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=4, max=10),
    retry=retry_if_exception_type(
        (requests.HTTPError, requests.ConnectionError, requests.Timeout)
    ),
)
def progress_retrieve(url, file, disable=False):
    headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}

    Path(file).parent.mkdir(parents=True, exist_ok=True)

    # Raise HTTPError for transient errors
    # 429: Too Many Requests (rate limiting)
    # 500, 502, 503, 504: Server errors
    response = requests.get(url, headers=headers, stream=True)
    if response.status_code in (429, 500, 502, 503, 504):
        response.raise_for_status()
    total_size = int(response.headers.get("content-length", 0))
    chunk_size = 1024

    with tqdm(
        total=total_size,
        unit="B",
        unit_scale=True,
        unit_divisor=1024,
        desc=str(file),
        disable=disable,
    ) as t:
        with open(file, "wb") as f:
            for data in response.iter_content(chunk_size=chunk_size):
                f.write(data)
                t.update(len(data))


def retry(func: Callable) -> Callable:
    """
    Retry decorator to run retry function on specific exceptions, before raising them.

    Can for example be used for debugging issues which are hard to replicate or
    for for handling retrieval errors.

    Currently catches:
    - fiona.errors.DriverError

    Parameters
    ----------
    retries : int
        Number of retries before raising the exception.
    delay : int
        Delay between retries in seconds.

    Returns
    -------
    callable
        A decorator function that can be used to wrap the function to be retried.
    """
    retries = 3
    delay = 5

    @wraps(func)
    def wrapper(*args, **kwargs):
        for attempt in range(retries):
            try:
                return func(*args, **kwargs)
            except fiona.errors.DriverError as e:
                logger.warning(
                    f"Attempt {attempt + 1} failed: {type(e).__name__} - {e}. "
                    f"Retrying..."
                )
                time.sleep(delay)
        raise Exception("Retrieval retries exhausted.")

    return wrapper


def mock_snakemake(
    rulename,
    root_dir=None,
    configfiles=None,
    submodule_dir="workflow/submodules/pypsa-eur",
    **wildcards,
):
    """
    This function is expected to be executed from the 'scripts'-directory of '
    the snakemake project. It returns a snakemake.script.Snakemake object,
    based on the Snakefile.

    If a rule has wildcards, you have to specify them in **wildcards.

    Parameters
    ----------
    rulename: str
        name of the rule for which the snakemake object should be generated
    root_dir: str/path-like
        path to the root directory of the snakemake project
    configfiles: list, str
        list of configfiles to be used to update the config
    submodule_dir: str, Path
        in case PyPSA-Eur is used as a submodule, submodule_dir is
        the path of pypsa-eur relative to the project directory.
    **wildcards:
        keyword arguments fixing the wildcards. Only necessary if wildcards are
        needed.
    """
    import os

    import snakemake as sm
    from pypsa.definitions.structures import Dict
    from snakemake.api import Workflow
    from snakemake.common import SNAKEFILE_CHOICES
    from snakemake.script import Snakemake
    from snakemake.settings.types import (
        ConfigSettings,
        DAGSettings,
        ResourceSettings,
        StorageSettings,
        WorkflowSettings,
    )

    script_dir = Path(__file__).parent.resolve()
    if root_dir is None:
        root_dir = script_dir.parent
    else:
        root_dir = Path(root_dir).resolve()

    workdir = None
    user_in_script_dir = Path.cwd().resolve() == script_dir
    if str(submodule_dir) in __file__:
        # the submodule_dir path is only need to locate the project dir
        os.chdir(Path(__file__[: __file__.find(str(submodule_dir))]))
    elif user_in_script_dir:
        os.chdir(root_dir)
    elif Path.cwd().resolve() != root_dir:
        logger.info(
            "Not in scripts or root directory, will assume this is a separate workdir"
        )
        workdir = Path.cwd()

    try:
        for p in SNAKEFILE_CHOICES:
            p = root_dir / p
            if os.path.exists(p):
                snakefile = p
                break
        if configfiles is None:
            configfiles = []
        elif isinstance(configfiles, str):
            configfiles = [configfiles]

        resource_settings = ResourceSettings()
        config_settings = ConfigSettings(configfiles=map(Path, configfiles))
        workflow_settings = WorkflowSettings()
        storage_settings = StorageSettings()
        dag_settings = DAGSettings(rerun_triggers=[])
        workflow = Workflow(
            config_settings,
            resource_settings,
            workflow_settings,
            storage_settings,
            dag_settings,
            storage_provider_settings=dict(),
            overwrite_workdir=workdir,
        )
        workflow.include(snakefile)

        if configfiles:
            for f in configfiles:
                if not os.path.exists(f):
                    raise FileNotFoundError(f"Config file {f} does not exist.")
                workflow.configfile(f)

        workflow.global_resources = {}
        rule = workflow.get_rule(rulename)
        dag = sm.dag.DAG(workflow, rules=[rule])
        wc = Dict(wildcards)
        job = sm.jobs.Job(rule, dag, wc)

        def make_accessable(*ios):
            for io in ios:
                for i, _ in enumerate(io):
                    io[i] = os.path.abspath(io[i])

        make_accessable(job.input, job.output, job.log)
        snakemake = Snakemake(
            job.input,
            job.output,
            job.params,
            job.wildcards,
            job.threads,
            job.resources,
            job.log,
            job.dag.workflow.config,
            job.rule.name,
            None,
        )
        # create log and output dir if not existent
        for path in list(snakemake.log) + list(snakemake.output):
            Path(path).parent.mkdir(parents=True, exist_ok=True)

    finally:
        if user_in_script_dir:
            os.chdir(script_dir)
    return snakemake


def generate_periodic_profiles(dt_index, nodes, weekly_profile, localize=None):
    """
    Give a 24*7 long list of weekly hourly profiles, generate this for each
    country for the period dt_index, taking account of time zones and summer
    time.
    """
    weekly_profile = pd.Series(weekly_profile, range(24 * 7))

    week_df = pd.DataFrame(index=dt_index, columns=nodes)

    for node in nodes:
        ct = node[:2] if node[:2] != "XK" else "RS"
        timezone = pytz.timezone(pytz.country_timezones[ct][0])
        tz_dt_index = dt_index.tz_convert(timezone)
        week_df[node] = [24 * dt.weekday() + dt.hour for dt in tz_dt_index]
        week_df[node] = week_df[node].map(weekly_profile)

    week_df = week_df.tz_localize(localize)

    return week_df


def parse(infix):
    """
    Recursively parse a chained wildcard expression into a dictionary or a YAML
    object.

    Parameters
    ----------
    list_to_parse : list
        The list to parse.

    Returns
    -------
    dict or YAML object
        The parsed list.
    """
    if len(infix) == 1:
        return yaml.safe_load(infix[0])
    else:
        return {infix.pop(0): parse(infix)}


def update_config_from_wildcards(config, w, inplace=True):
    """
    Parses configuration settings from wildcards and updates the config.
    """

    if not inplace:
        config = copy.deepcopy(config)

    if w.get("opts"):
        opts = w.opts.split("-")

        if nhours := get_opt(opts, r"^\d+(h|seg)$"):
            config["clustering"]["temporal"]["resolution_elec"] = nhours

        co2l_enable, co2l_value = find_opt(opts, "Co2L")
        if co2l_enable:
            config["electricity"]["co2limit_enable"] = True
            if co2l_value is not None:
                config["electricity"]["co2limit"] = (
                    co2l_value * config["electricity"]["co2base"]
                )

        gasl_enable, gasl_value = find_opt(opts, "CH4L")
        if gasl_enable:
            config["electricity"]["gaslimit_enable"] = True
            if gasl_value is not None:
                config["electricity"]["gaslimit"] = gasl_value * 1e6

        if "Ept" in opts:
            config["costs"]["emission_prices"]["co2_monthly_prices"] = True

        ep_enable, ep_value = find_opt(opts, "Ep")
        if ep_enable:
            config["costs"]["emission_prices"]["enable"] = True
            if ep_value is not None:
                config["costs"]["emission_prices"]["co2"] = ep_value

        if "ATK" in opts:
            config["autarky"]["enable"] = True
            if "ATKc" in opts:
                config["autarky"]["by_country"] = True

        attr_lookup = {
            "p": "p_nom_max",
            "e": "e_nom_max",
            "c": "capital_cost",
            "m": "marginal_cost",
        }
        for o in opts:
            flags = ["+e", "+p", "+m", "+c"]
            if all(flag not in o for flag in flags):
                continue
            carrier, component, attr_factor = o.split("+")
            attr = attr_lookup[attr_factor[0]]
            factor = float(attr_factor[1:])
            if not isinstance(config["adjustments"]["electricity"], dict):
                config["adjustments"]["electricity"] = dict()
            update_config(
                config["adjustments"]["electricity"],
                {"factor": {component: {carrier: {attr: factor}}}},
            )

        for o in opts:
            if o.startswith("lv") or o.startswith("lc"):
                config["electricity"]["transmission_limit"] = o[1:]
                break

    if w.get("sector_opts"):
        opts = w.sector_opts.split("-")

        if "T" in opts:
            config["sector"]["transport"] = True

        if "H" in opts:
            config["sector"]["heating"] = True

        if "B" in opts:
            config["sector"]["biomass"] = True

        if "I" in opts:
            config["sector"]["industry"] = True

        if "A" in opts:
            config["sector"]["agriculture"] = True

        if "CCL" in opts:
            config["solving"]["constraints"]["CCL"] = True

        eq_value = get_opt(opts, r"^EQ+\d*\.?\d+(c|)")
        for o in opts:
            if eq_value is not None:
                config["solving"]["constraints"]["EQ"] = eq_value
            elif "EQ" in o:
                config["solving"]["constraints"]["EQ"] = True
            break

        if "BAU" in opts:
            config["solving"]["constraints"]["BAU"] = True

        if "SAFE" in opts:
            config["solving"]["constraints"]["SAFE"] = True

        if nhours := get_opt(opts, r"^\d+(h|sn|seg)$"):
            config["clustering"]["temporal"]["resolution_sector"] = nhours

        if "decentral" in opts:
            config["sector"]["electricity_transmission_grid"] = False

        if "noH2network" in opts:
            config["sector"]["H2_network"] = False

        if "nowasteheat" in opts:
            config["sector"]["use_fischer_tropsch_waste_heat"] = False
            config["sector"]["use_methanolisation_waste_heat"] = False
            config["sector"]["use_haber_bosch_waste_heat"] = False
            config["sector"]["use_methanation_waste_heat"] = False
            config["sector"]["use_fuel_cell_waste_heat"] = False
            config["sector"]["use_electrolysis_waste_heat"] = False

        if "nodistrict" in opts:
            config["sector"]["district_heating"]["progress"] = 0.0

        dg_enable, dg_factor = find_opt(opts, "dist")
        if dg_enable:
            config["sector"]["electricity_distribution_grid"] = True
            if dg_factor is not None:
                config["sector"]["electricity_distribution_grid_cost_factor"] = (
                    dg_factor
                )

        if "biomasstransport" in opts:
            config["sector"]["biomass_transport"] = True

        _, maxext = find_opt(opts, "linemaxext")
        if maxext is not None:
            config["lines"]["max_extension"] = maxext * 1e3
            config["links"]["max_extension"] = maxext * 1e3

        _, co2l_value = find_opt(opts, "Co2L")
        if co2l_value is not None:
            config["co2_budget"] = float(co2l_value)

        if co2_distribution := get_opt(opts, r"^(cb)\d+(\.\d+)?(ex|be)$"):
            config["co2_budget"] = co2_distribution

        if co2_budget := get_opt(opts, r"^(cb)\d+(\.\d+)?$"):
            config["co2_budget"] = float(co2_budget[2:])

        attr_lookup = {
            "p": "p_nom_max",
            "e": "e_nom_max",
            "c": "capital_cost",
            "m": "marginal_cost",
        }
        for o in opts:
            flags = ["+e", "+p", "+m", "+c"]
            if all(flag not in o for flag in flags):
                continue
            carrier, component, attr_factor = o.split("+")
            attr = attr_lookup[attr_factor[0]]
            factor = float(attr_factor[1:])
            if not isinstance(config["adjustments"]["sector"], dict):
                config["adjustments"]["sector"] = dict()
            update_config(
                config["adjustments"]["sector"],
                {"factor": {component: {carrier: {attr: factor}}}},
            )

        _, sdr_value = find_opt(opts, "sdr")
        if sdr_value is not None:
            config["costs"]["social_discountrate"] = sdr_value / 100

        _, seq_limit = find_opt(opts, "seq")
        if seq_limit is not None:
            config["sector"]["co2_sequestration_potential"] = seq_limit

        # any config option can be represented in wildcard
        for o in opts:
            if o.startswith("CF+"):
                infix = o.split("+")[1:]
                update_config(config, parse(infix))

    if not inplace:
        return config


@tenacity_retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=4, max=10),
    retry=retry_if_exception_type(
        (requests.HTTPError, requests.ConnectionError, requests.Timeout)
    ),
)
def get_checksum_from_zenodo(file_url):
    parts = file_url.split("/")
    record_id = parts[parts.index("records") + 1]
    filename = parts[-1]

    response = requests.get(f"https://zenodo.org/api/records/{record_id}", timeout=30)
    # Raise HTTPError for transient errors
    # 429: Too Many Requests (rate limiting)
    # 500, 502, 503, 504: Server errors
    if response.status_code in (429, 500, 502, 503, 504):
        response.raise_for_status()
    response.raise_for_status()
    data = response.json()

    for file in data["files"]:
        if file["key"] == filename:
            return file["checksum"]
    return None


def validate_checksum(file_path, zenodo_url=None, checksum=None):
    """
    Validate file checksum against provided or Zenodo-retrieved checksum.
    Calculates the hash of a file using 64KB chunks. Compares it against a
    given checksum or one from a Zenodo URL.

    Parameters
    ----------
    file_path : str
        Path to the file for checksum validation.
    zenodo_url : str, optional
        URL of the file on Zenodo to fetch the checksum.
    checksum : str, optional
        Checksum (format 'hash_type:checksum_value') for validation.

    Raises
    ------
    AssertionError
        If the checksum does not match, or if neither `checksum` nor `zenodo_url` is provided.


    Examples
    --------
    >>> validate_checksum("/path/to/file", checksum="md5:abc123...")
    >>> validate_checksum(
    ...     "/path/to/file",
    ...     zenodo_url="https://zenodo.org/records/12345/files/example.txt",
    ... )

    If the checksum is invalid, an AssertionError will be raised.
    """
    assert checksum or zenodo_url, "Either checksum or zenodo_url must be provided"
    if zenodo_url:
        checksum = get_checksum_from_zenodo(zenodo_url)
    hash_type, checksum = checksum.split(":")
    hasher = hashlib.new(hash_type)
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):  # 64kb chunks
            hasher.update(chunk)
    calculated_checksum = hasher.hexdigest()
    assert calculated_checksum == checksum, (
        "Checksum is invalid. This may be due to an incomplete download. Delete the file and re-execute the rule."
    )


def get_snapshots(
    snapshots: dict, drop_leap_day: bool = False, freq: str = "h", **kwargs
) -> pd.DatetimeIndex:
    """
    Returns a DateTimeIndex of snapshots, supporting multiple time ranges.

    Parameters
    ----------
    snapshots : dict
        Dictionary containing time range parameters. 'start' and 'end' can be
        strings or lists of strings for multiple date ranges.
    drop_leap_day : bool, default False
        If True, removes February 29th from the DateTimeIndex in leap years.
    freq : str, default "h"
        Frequency string indicating the time step interval (e.g., "h" for hourly)
    **kwargs : dict
        Additional keyword arguments passed to pd.date_range().

    Returns
    -------
    pd.DatetimeIndex
    """
    start = (
        snapshots["start"]
        if isinstance(snapshots["start"], list)
        else [snapshots["start"]]
    )
    end = snapshots["end"] if isinstance(snapshots["end"], list) else [snapshots["end"]]

    assert len(start) == len(end), (
        "Lists of start and end dates must have the same length"
    )

    time_periods = []
    for s, e in zip(start, end):
        period = pd.date_range(
            start=s, end=e, freq=freq, inclusive=snapshots["inclusive"], **kwargs
        )
        time_periods.append(period)

    time = pd.DatetimeIndex([])
    for period in time_periods:
        time = time.append(period)

    if drop_leap_day and time.is_leap_year.any():
        time = time[~((time.month == 2) & (time.day == 29))]

    return time


def sanitize_custom_columns(n: pypsa.Network):
    """
    Sanitize non-standard columns used throughout the workflow.

    Parameters
    ----------
        n (pypsa.Network): The network object.

    Returns
    -------
        None
    """
    if "reversed" in n.links.columns:
        # Replace NA values with default value False
        n.links.loc[n.links.reversed.isna(), "reversed"] = False
        n.links.reversed = n.links.reversed.astype(bool)


def rename_techs(label: str) -> str:
    """
    Rename technology labels for better readability.

    Removes some prefixes and renames if certain conditions defined in function body are met.

    Parameters
    ----------
    label: str
        Technology label to be renamed

    Returns
    -------
    str
        Renamed label
    """
    prefix_to_remove = [
        "residential ",
        "services ",
        "urban ",
        "rural ",
        "central ",
        "decentral ",
    ]

    rename_if_contains = [
        "CHP",
        "gas boiler",
        "biogas",
        "solar thermal",
        "air heat pump",
        "ground heat pump",
        "resistive heater",
        "Fischer-Tropsch",
    ]

    rename_if_contains_dict = {
        "water tanks": "hot water storage",
        "retrofitting": "building retrofitting",
        # "H2 Electrolysis": "hydrogen storage",
        # "H2 Fuel Cell": "hydrogen storage",
        # "H2 pipeline": "hydrogen storage",
        "battery": "battery storage",
        "H2 for industry": "H2 for industry",
        "land transport fuel cell": "land transport fuel cell",
        "land transport oil": "land transport oil",
        "oil shipping": "shipping oil",
        # "CC": "CC"
    }

    rename = {
        "solar": "solar PV",
        "Sabatier": "methanation",
        "offwind": "offshore wind",
        "offwind-ac": "offshore wind (AC)",
        "offwind-dc": "offshore wind (DC)",
        "offwind-float": "offshore wind (Float)",
        "onwind": "onshore wind",
        "ror": "hydroelectricity",
        "hydro": "hydroelectricity",
        "PHS": "hydroelectricity",
        "NH3": "ammonia",
        "co2 Store": "DAC",
        "co2 stored": "CO2 sequestration",
        "AC": "transmission lines",
        "DC": "transmission lines",
        "B2B": "transmission lines",
    }

    for ptr in prefix_to_remove:
        if label[: len(ptr)] == ptr:
            label = label[len(ptr) :]

    for rif in rename_if_contains:
        if rif in label:
            label = rif

    for old, new in rename_if_contains_dict.items():
        if old in label:
            label = new

    for old, new in rename.items():
        if old == label:
            label = new
    return label


def load_cutout(
    cutout_files: str | list[str], time: None | pd.DatetimeIndex = None
) -> atlite.Cutout:
    """
    Load and optionally combine multiple cutout files.

    Parameters
    ----------
    cutout_files : str or list of str
        Path to a single cutout file or a list of paths to multiple cutout files.
        If a list is provided, the cutouts will be concatenated along the time dimension.
    time : pd.DatetimeIndex, optional
        If provided, select only the specified times from the cutout.

    Returns
    -------
    atlite.Cutout
        Merged cutout with optional time selection applied.
    """
    if isinstance(cutout_files, str):
        cutout = atlite.Cutout(cutout_files)
    elif isinstance(cutout_files, list):
        cutout_da = [atlite.Cutout(c).data for c in cutout_files]
        combined_data = xr.concat(cutout_da, dim="time", data_vars="minimal")
        cutout = atlite.Cutout(NamedTemporaryFile().name, data=combined_data)

    if time is not None:
        cutout.data = cutout.data.sel(time=time)

    return cutout


def load_costs(cost_file: str) -> pd.DataFrame:
    """
    Load prepared cost data from CSV.

    Parameters
    ----------
    cost_file : str
        Path to the CSV file containing cost data

    Returns
    -------
    costs : pd.DataFrame
        DataFrame containing the prepared cost data
    """

    return pd.read_csv(cost_file, index_col=0)


def make_index(
    c, cname0="bus0", cname1="bus1", prefix="", connector="->", suffix="", separator=" "
):
    idx = [prefix, c[cname0], connector, c[cname1], suffix]
    idx = [i for i in idx if i]
    return separator.join(idx)


def extract_grid_data_tyndp(
    links,
    replace_dict: dict = {},
    expand_from_index: bool = True,
    idx_prefix: str = "",
    idx_connector: str = "",
    idx_suffix: str = "",
    idx_separator: str = " ",
):
    """
    Extract TYNDP reference grid data from the raw input table.

    Parameters
    ----------
    links : pd.DataFrame
        DataFrame with raw links to extract grid information from
    replace_dict : dict
        Dictionary with region names to replace
    expand_from_index : bool
        Whether to expand the bus0 and bus1 from index or directly use the columns
    idx_prefix : str, optional
        Prefix to prepend to generated indices.
    idx_connector : str, optional
        Separator string between bus0 and bus1 in generated indices (e.g., "->", "-").
    idx_suffix : str, optional
        Suffix to append to generated indices.
    idx_separator : str, optional
        String used to join index components (prefix, bus0, connector, bus1, suffix).

    Returns
    -------
    pd.DataFrame
        DataFrame with extracted grid data information with nominal capacity in input unit, bus0 and bus1
    """

    if expand_from_index:
        links.loc[:, "Border"] = links["Border"].replace(replace_dict, regex=True)
        links = pd.concat(
            [
                links,
                links.Border.str.split("-", expand=True).set_axis(
                    ["bus0", "bus1"], axis=1
                ),
            ],
            axis=1,
        )
    elif "bus0" in links.columns and "bus1" in links.columns:
        links.loc[:, ["bus0", "bus1"]] = links[["bus0", "bus1"]].replace(
            replace_dict, regex=True
        )
    else:
        raise KeyError(
            f"Columns 'bus0' and 'bus1' must be present in the input DataFrame. Available columns: {list(links.columns)}"
        )

    # Create forward and reverse direction dataframes
    forward_links = links[["bus0", "bus1", "Summary Direction 1"]].rename(
        columns={"Summary Direction 1": "p_nom"}
    )

    reverse_links = links[["bus1", "bus0", "Summary Direction 2"]].rename(
        columns={"bus1": "bus0", "bus0": "bus1", "Summary Direction 2": "p_nom"}
    )

    # Combine into unidirectional links and return
    links = pd.concat([forward_links, reverse_links])

    links.index = links.apply(
        make_index,
        axis=1,
        prefix=idx_prefix,
        connector=idx_connector,
        suffix=idx_suffix,
        separator=idx_separator,
    )

    return links


def safe_pyear(
    year: int | str,
    available_years: list = [2030, 2040, 2050],
    source: str = "TYNDP",
    verbose: bool = True,
) -> int:
    """
    Checks and adjusts whether a given pyear is in the available years of a given data source. If not, it
    falls back to the previous available year.

    Parameters
    ----------
    year : int
        Planning horizon year which will be checked and possibly adjusted to previous available year.
    available_years : list, optional
        List of available years. Defaults to [2030, 2040, 2050].
    source : str, optional
        Source of the data for which availability will be checked. For logging purpose only. Defaults to "TYNDP".
    verbose : bool, optional
        Whether to activate verbose logging. Defaults to True.

    Returns
    -------
    year_new : int
        Safe pyear adjusted for available years
    """

    if not available_years:
        raise ValueError(
            "No `available_years` provided. Expected a non-empty list of years."
        )
    if not isinstance(year, int):
        year = int(year)
    if year not in available_years:
        year_new = available_years[
            bisect_right(sorted(available_years), year, lo=1) - 1
        ]
        if verbose:
            logger.warning(
                f"{source} data unavailable for planning horizon {year}. Falling back to previous available year {year_new}."
            )
    else:
        year_new = year

    return year_new


def map_tyndp_carrier_names(
    df: pd.DataFrame,
    carrier_mapping_fn: str,
    on_columns: list[str],
    drop_on_columns=False,
):
    """
    Map external carriers to available tyndp_carrier names based on an input mapping. Optionally drop merged on columns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with external carriers to map
    carrier_mapping_fn : str
        Path to file with mapping from external carriers to available tyndp_carrier names.
    on_columns : list[str]
        Columns to merge on between the external carriers and tyndp_carriers.
    drop_on_columns : bool, optional
        Whether to drop merge columns and rename `open_tyndp_carrier` and `open_tyndp_index` to `carrier`
        and `index_carrier`. Defaults to False.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with external carriers mapped to available tyndp_carriers and index_carriers.
    """

    # Read TYNDP carrier mapping
    carrier_mapping = (
        pd.read_csv(carrier_mapping_fn)[
            on_columns + ["open_tyndp_carrier", "open_tyndp_index", "open_tyndp_type"]
        ]
    ).dropna()

    # Map the carriers
    df = df.merge(carrier_mapping, on=on_columns, how="left")

    # If the carrier is DSR or Other Non-RES, the different price bands are too diverse for a robust external
    # mapping. Instead, we will combine the carrier and type information.
    if "pemmdb_carrier" in on_columns:

        def normalize_carrier(s):
            return s.lower().replace(" ", "-").replace("other-non-res", "chp")

        # Other Non-RES are assumed to represent CHP plants (according to TYNDP 2024 Methodology report p.37)
        df = df.assign(
            open_tyndp_carrier=lambda x: np.where(
                x["pemmdb_carrier"].isin(["DSR", "Other Non-RES"]),
                x["pemmdb_carrier"].apply(normalize_carrier),
                x["open_tyndp_carrier"],
            ),
            open_tyndp_index=lambda x: np.where(
                x["pemmdb_carrier"].isin(["DSR", "Other Non-RES"]),
                x["open_tyndp_carrier"]
                + "-"
                + x["pemmdb_type"].apply(normalize_carrier),
                x["open_tyndp_index"],
            ),
        )

    if not drop_on_columns:
        return df

    # Otherwise drop merge columns and rename to new "carrier" and "index_carrier" column
    df = df.drop(on_columns, axis="columns").rename(
        columns={
            "open_tyndp_carrier": "carrier",
            "open_tyndp_index": "index_carrier",
        }
    )

    # Move "carrier" and "index_carrier" to the front
    cols = ["carrier", "index_carrier"] + [
        col for col in df.columns if col not in ["carrier", "index_carrier"]
    ]

    return df[cols]


def get_version(hash_len: int = 9) -> str:
    """
    Create a version identifier from git repository state.

    Returns a version string based on the latest reachable tag and current commit:
    - If HEAD is exactly at a tag: returns the tag name (e.g., "v1.2.3")
    - If HEAD is beyond a tag: returns "tag+g{hash}" (e.g., "v1.2.3+g1a2b3c4d")
    - If no tags found: returns just the commit hash (e.g., "1a2b3c4d5")
    """
    try:
        repo = git.Repo(search_parent_directories=True)
        tags = sorted(
            repo.tags, key=lambda t: t.commit.committed_datetime, reverse=True
        )
        last_tag = None
        for tag in tags:
            if repo.is_ancestor(tag.commit, repo.head.commit):
                last_tag = tag
                break
        if last_tag and last_tag.commit == repo.head.commit:
            return f"{last_tag}"
        elif last_tag:
            return f"{last_tag}+g{repo.head.commit.hexsha[:hash_len]}"
        else:
            return repo.head.commit.hexsha[:hash_len]

    except Exception as e:
        logger.warning(f"Failed to determine version from git repository: {e}")
        return "unknown"


def convert_units(
    df: pd.DataFrame,
    unit_col: str = "unit",
    value_col: str = "value",
    invert: bool = False,
) -> pd.DataFrame:
    """
    Convert values to standardized units based on unit type.

    Energy values are converted to MWh, power values to MW:
    - Energy units (TWh, GWh, MWh, kWh) → MWh
    - Power units (GW, MW, kW) → MW

    When invert=False (default):
        - Values are converted from unit_col units to standard units (MWh/MW)
        - The "unit" column is updated to reflect the standardized unit

    When invert=True:
        - Values are converted from standard units (MWh/MW) back to unit_col units
        - The "unit" column is NOT modified
        - Useful for reverting previously standardized data

    Parameters
    ----------
    df : pd.DataFrame
        Long-format DataFrame containing values to convert.
    unit_col : str, default "unit"
        Name of the column containing the unit information.
        When invert=False: contains source units to convert from.
        When invert=True: contains target units to convert to.
    value_col : str, default "value"
        Name of the column containing values to convert.
    invert : bool, default False
        If False, convert to standard units and update "unit" column.
        If True, convert from standard units using inverse factors without modifying "unit" column.

    Returns
    -------
    pd.DataFrame
        DataFrame with converted values.
    """
    df = df.copy()

    unit_conversion = {
        "TWh": 1000000,
        "GWh": 1000,
        "MWh": 1,
        "GW": 1000,
        "MW": 1,
        "kW": 0.001,
    }

    if invert:
        # Inverse conversion factor to revert unit
        unit_conversion = {k: 1 / v for k, v in unit_conversion.items()}

    # Convert values using conversion factors
    conversion_factors = df[unit_col].map(unit_conversion)
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce") * conversion_factors

    # Update unit column
    if not invert:
        df["unit"] = df[unit_col].apply(
            lambda x: "MWh" if x in ENERGY_UNITS else "MW" if x in POWER_UNITS else x
        )

    return df


def check_cyear(cyear: int, scenario: str) -> int:
    """Check if the climatic year is valid for the given scenario."""

    valid_years = {
        "NT": np.arange(1983, 2018).tolist(),
        "DE": [1995, 2008, 2009],
        "GA": [1995, 2008, 2009],
    }

    if cyear not in valid_years[scenario]:
        logger.warning(
            f"Snapshot year {cyear} doesn't match available TYNDP data. Falling back to 2009."
        )
        cyear = 2009

    return cyear


def get_tyndp_conventional_thermals(
    mapping: pd.DataFrame,
    tyndp_conventional_carriers: list[str],
    group_conventionals: bool,
    include_h2_fuel_cell: bool,
    include_h2_turbine: bool,
) -> tuple[dict[str, str], list[str]]:
    """
    Get list of TYNDP conventional thermal generation technologies.

    Parameters
    ----------
    mapping : pd.DataFrame
        TYNDP conventional carrier mapping (grouped or ungrouped).
    tyndp_conventional_carriers : list[str]
        TYNDP conventional carriers.
    group_conventionals : bool
        Whether to group conventional thermal technologies.
    include_h2_fuel_cell : bool
        Whether to include hydrogen fuel cell technology.
    include_h2_turbine : bool
        Whether to include hydrogen turbine technology.

    Returns
    -------
    tuple[dict[str, str], list[str]]
        Dictionary with conventional mapping and list of conventional thermal technology names.
    """

    # Filter mapping for conventional carriers while setting oil as common fuel for oil technologies
    mapping = (
        mapping[["open_tyndp_carrier", "open_tyndp_type", "pypsa_eur_carrier"]]
        .query("open_tyndp_carrier in @tyndp_conventional_carriers")
        .replace(
            {"open_tyndp_carrier": ["oil-light", "oil-heavy", "oil-shale"]}, "oil"
        )  # TODO To remove once the three carriers have been implemented
    )

    if group_conventionals:
        mapping = mapping.groupby("open_tyndp_type").first()

    conventional_dict = mapping.open_tyndp_carrier.to_dict()
    conventional_thermals = list(conventional_dict)

    if include_h2_fuel_cell:
        conventional_thermals.append("h2-fuel-cell")
    if include_h2_turbine:
        conventional_thermals.append("h2-ccgt")

    return conventional_dict, conventional_thermals


def interpolate_demand(
    available_years: list[int],
    pyear: int,
    load_single_year_func: Callable,
    **load_kwargs,
) -> pd.DataFrame | pd.Series:
    """
    Interpolate demand between available years.

    Parameters
    ----------
    available_years : list[int]
        Sorted list of years for which data is available.
    pyear : int
        Planning year to interpolate demand for.
    load_single_year_func : Callable
        Function to load data for a single planning year.
    **load_kwargs
        Keyword arguments to pass to load_single_year_func. Must include 'pyear'
        as a parameter key, which will be overridden with interpolation boundary years.

    Returns
    -------
    pd.DataFrame | pd.Series
        Interpolated demand data.
    """
    # Currently, only interpolation is implemented, not extrapolation
    idx = bisect_right(available_years, pyear)
    if idx == 0:
        # Planning horizon is before all available years
        logger.warning(
            f"Year {pyear} is before the first available year {available_years[0]}. "
            f"Falling back to first available year."
        )
        year_lower = year_upper = available_years[0]
    elif idx == len(available_years):
        # Planning horizon is after all available years
        logger.warning(
            f"Year {pyear} is after the latest available year {available_years[-1]}. "
            f"Falling back to latest available year."
        )
        year_lower = year_upper = available_years[-1]
    else:
        year_lower = available_years[idx - 1]
        year_upper = available_years[idx]

    logger.debug(f"Interpolating {pyear} from {year_lower} and {year_upper}")

    kwargs_lower = {**load_kwargs, "pyear": year_lower}
    kwargs_upper = {**load_kwargs, "pyear": year_upper}

    df_lower = load_single_year_func(**kwargs_lower)
    df_upper = load_single_year_func(**kwargs_upper)

    # Check if data was loaded successfully
    if df_lower.empty and df_upper.empty:
        logger.error("Both years failed to load")
        return pd.DataFrame()
    elif df_lower.empty:
        logger.warning(
            f"Year {year_lower} failed to load. Filling with zeros for interpolation."
        )
        df_lower = pd.DataFrame(0, index=df_upper.index, columns=df_upper.columns)
    elif df_upper.empty:
        logger.warning(
            f"Year {year_upper} failed to load. Using data from lower year for interpolation."
        )
        df_upper = df_lower

    if year_upper == year_lower:
        return df_lower

    # Handle column mismatches for DataFrames (only relevant for DataFrame, not Series)
    if isinstance(df_lower, pd.DataFrame) and isinstance(df_upper, pd.DataFrame):
        missing_in_lower = df_upper.columns.difference(df_lower.columns)
        missing_in_upper = df_lower.columns.difference(df_upper.columns)

        if len(missing_in_lower) > 0 or len(missing_in_upper) > 0:
            logger.warning(
                f"Column mismatch between {year_lower} and {year_upper}. "
                f"Missing columns filled with zeros. "
                f"Missing in {year_lower}: {list(missing_in_lower)}, "
                f"Missing in {year_upper}: {list(missing_in_upper)}"
            )
        df_lower_aligned, df_upper_aligned = df_lower.align(
            df_upper, join="outer", axis=1, fill_value=0
        )
    else:
        # For Series, just align
        df_lower_aligned, df_upper_aligned = df_lower.align(
            df_upper, join="outer", fill_value=0
        )

    # Perform linear interpolation
    weight = (pyear - year_lower) / (year_upper - year_lower)
    result = df_lower_aligned * (1 - weight) + df_upper_aligned * weight

    return result
