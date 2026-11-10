<!-- SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp> -->
<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!---->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Retrieving Data {#data}

Not all data dependencies are shipped with the git repository, since git is not suited for handling large changing files.
Instead we use separate steps in the workflow (`rules` executed by `snakemake`) to download external data using the `retrieve_<dataset>` rules.

Data is generally retrieved in a version-controlled manner, enabling control over input data versions, reproducibility and consistency of modelling runs.
The rules download data into subfolders in the `data/` directory, following the structure
`data/{dataset}/{source}/{version}`, e.g. `data/jrc_idees/primary/March-2025-V1/`.
Which specific data version is retrieved can be controlled in the [data configuration](configuration.md#data_cf).

Every dataset is listed with its owner, link and license in the
[data inventory](data_sources.md#data-inventory). Its available versions and
sources (`archive`, `primary` or `build`) are registered in `data/versions.csv`.
Most `retrieve_<dataset>` rules simply download the registered URL for the
selected source. The rules below need more than that.

For Open-TYNDP runs, most datasets can also be retrieved from a dedicated Google Cloud Storage
bucket instead of their original sources. See [tyndp_archive](sb.md#tyndp_archive) in the SB documentation.

## Local data cache {#local_cache}

The local cache keeps every retrieved file in a single source-agnostic tree, which makes it
possible to run the workflow on a machine with no internet access at all. It drops the
`{source}` segment from the paths above, so datasets are held once regardless of where it was
fetched from under `{local_cache.directory}/{dataset}/{version}`.

Setting [`data: local_cache: enable:`](configuration.md#data_cf) points every dataset at that
tree and disables retrieval. Per default the local cache directory is set as `data/local-cache`
but can also take on any other relative or absolute path the user sets. 
The relevant configuration block reads:

```yaml
data:
  local_cache:
    enable: true
    directory: data/local-cache
    fill: false
```

### Filling the cache with network access {#local_cache_online}

The cache is filled in a separate step, depending on the workflow. The CBA workflow rule will, however,
also retrieve all required files to run the prerequisite SB workflow, unless the SB networks are taken pre-solved (see below).

```console
# everything Scenario Building needs
$ python utils/collect_data.py --configfile path/to/your/config.yaml
# or
$ pixi run collect-data

# everything the CBA needs + required SB retrieves
$ python utils/collect_data.py --cba --configfile path/to/your/config.yaml
# or
$ pixi run collect-data-cba
```

Both tasks run `utils/collect_data.py`. The script dry-runs the full Scenario Building graph,
with `--forceall`, so every `retrieve_*` rule that graph contains is listed based on your config settings
regardless of what a previous run already left in `resources/`. Then everything missing from the cache
is downloaded. Anything else passed on the command line goes through to Snakemake, so
`pixi run collect-data -n` lists what would be downloaded without fetching anything and
`pixi run collect-data --configfile config/test/config.tyndp.yaml` collects for another config.
The dry runs that only work out which datasets are needed stay quiet unless they fail; add the
`--verbose` flag to see them, and every other Snakemake call, in full.
The source selection still applies: with `data_config: tyndp` the cache fills from the
Google Cloud Storage mirror, into the same paths.

!!! note "Collect data without pixi"
    `collect-data` and `collect-data-cba` are pixi tasks, not Snakemake rules — they just run the
    collect script via the commands below. If you don't use pixi, run them directly in any
    environment that has Snakemake and the project's other dependencies installed.

    Equivalent to `pixi run collect-data`:
    `python utils/collect_data.py --configfile config/config.tyndp.yaml`

    Equivalent to `pixi run collect-data-cba`:
    `python utils/collect_data.py --cba --configfile config/config.tyndp.yaml`

!!! note "How the CBA task handles the `clean_projects` checkpoint"
    Most of the CBA graph, the per-project networks and the Scenario Building chain they build
    on, only appears once the `clean_projects` checkpoint has run, and forcing the graph would
    keep it from running. `pixi run collect-data-cba` therefore starts with the
    [`collect_cba_data`](cba_rules.md#rule-collect_cba_data) rule, which gathers what the checkpoint
    reads, along with the pre-solved SB networks if configured, and runs the checkpoint itself.
    It then collects the Scenario Building graph, which covers everything the expanded CBA
    graph goes on to need, and checks afterwards that nothing is left, stopping with the rule names
    if a dataset is uncollected.

To reach an offline machine, fill the cache where the network is available and copy the
directory across, e.g. using rsync:

```console
$ rsync -a data/local-cache/ offline-machine:~/open-tyndp/data/local-cache/
```

!!! note "Repository data"
    The cache holds retrieved data only, so the offline machine needs the repository as well.
    Clone it where there is network and copy the whole checkout across, or clone it on the
    offline machine if it can reach the Git remote. `data/` matters here: Some files like
    `data/cba/*.csv`, `data/tyndp_technology_map.csv`, `data/tyndp_versions.csv` and the
    transmission project CSVs ship directly with the repository and are never fetched, so a copy
    that skips `data/` will break the workflow.

With `cba: cba_scenario_input: use_presolved: true` the pre-solved SB networks are collected into
`results/` rather than into the cache. The [`collect_cba_data`](cba_rules.md#rule-collect_cba_data)
rule lists them as an input, so the collect task retrieves them for the planning horizons the CBA
covers; copy `results/` across alongside the cache. That run still needs part of the Scenario
Building data, for the shapes, bidding zones and costs the CBA builds itself, so SB dependencies are
still retrieved into the local cache.

### Reading the cache without network access {#local_cache_offline}

A cached run contacts no remote URL. The retrieve rules stay defined but take the cache
manifest `{local_cache.directory}/.collected` as their input instead of a URL, because Snakemake
otherwise probes remote inputs over the network while the graph is built; filling the cache also drops
Snakemake's provenance records for the cached files, so that switch does not mark every retrieve
job out of date. The manifest and the provenance reset are written only when a collect task
finishes, so an interrupted `collect-data` leaves stale records behind and the next run will
want to re-fetch.

```console
# run either workflow, reading the cache if enabled, no network required
$ pixi run tyndp-sb
$ pixi run tyndp-cba
```

Should the cache be missing a file on offline execution, the workflow stops and
names the absent files, also on dry runs. That check reads the manifest, so it covers only what
the cache was told to hold. A dataset that was never collected is caught a moment later instead:
the run aborts before the first job, naming the retrieve rules that would have run.

Data already retrieved into `data/{dataset}/{source}/{version}` is not reused, so the first
`collect-data` downloads it again unless you move those folders into the cache layout by hand.

## Rules requiring credentials {#credentials}

Credentials are read from environment variables. You can put them into a `.env`
file in the repository root, which is loaded automatically and ignored by git.

| Rule | Requirement |
|------|-------------|
| [build_cutout][] | A [Copernicus Climate Data Store](https://cds.climate.copernicus.eu) account with the `cdsapi` key set up as described [on their website](https://cds.climate.copernicus.eu/how-to-api). Only needed when `data: cutout: source: build`. |
| `retrieve_electricity_demand_entsoe` | An [ENTSO-E Transparency Platform](https://transparency.entsoe.eu) API token in `ENTSOE_API_TOKEN`. Only needed when `data: entsoe_electricity_demand: source: build`. |
| `retrieve_corine` | A [Copernicus Land Monitoring Service](https://land.copernicus.eu/user/login) API key in `CORINE_API_TOKEN`. Only needed when `data: corine: source: primary`. |
| `retrieve_seawater_temperature` | [Copernicus Marine Service](https://marine.copernicus.eu/) credentials configured for the `copernicusmarine` package. Only needed when heat pumps use the `sea_water` heat source with a non-test cutout. |

## Rules with special handling {#special}

- `retrieve_cutout` downloads pre-built weather cutouts; see [cutouts](configuration.md#atlite_cf).
- `retrieve_osm_data_raw` queries the [Overpass API](https://overpass-api.de) per country when `data: osm: source: build`. The endpoint, retries and user agent are set under [overpass_api](configuration.md#overpass_api_cf).
- `retrieve_wdpa` and `retrieve_wdpa_marine` resolve the monthly changing download URL of the [World Database on Protected Areas](https://www.protectedplanet.net/). The data may not be redistributed, so the `archive` source points to a web archive copy.
- `retrieve_bidding_zones_entsoepy` and `retrieve_bidding_zones_electricitymaps` download bidding zone shapes via the [entsoe-py](https://github.com/EnergieID/entsoe-py) package and from [Electricity Maps](https://github.com/electricitymaps/electricitymaps-contrib). They are combined by [build_bidding_zones][].
- `retrieve_cost_data` downloads techno-economic assumptions from the [technology-data repository](https://github.com/pypsa/technology-data) as `data/costs/{source}/{version}/costs_{horizon}.csv`. The cost year can be fixed with `costs: year:` (see [costs](configuration.md#costs_cf)).

## Electricity demand data {#demand}

Historical hourly electricity demand is retrieved from three sources and
combined by [build_electricity_demand][]:

| Rule | Source | Output |
|------|--------|--------|
| `retrieve_electricity_demand_opsd` | [OPSD platform](https://data.open-power-system-data.org/time_series/2019-06-05/time_series_60min_singleindex.csv), per country | `data/opsd_electricity_demand/{source}/{version}/electricity_demand_opsd_raw.csv` |
| `retrieve_electricity_demand_entsoe` | [ENTSO-E Transparency Platform](https://transparency.entsoe.eu), per country | `data/entsoe_electricity_demand/{source}/{version}/electricity_demand_entsoe_raw.csv` |
| `retrieve_electricity_demand_neso` | [NESO Data Portal](https://www.neso.energy/data-portal/historic-demand-data), United Kingdom | `data/neso_electricity_demand/{source}/{version}/electricity_demand_neso_raw.csv` |

The spatial distribution of demand within countries uses:

- `retrieve_electricity_demand_energy_atlas`: a 1 km by 1 km raster of estimated annual electricity demand from the [JRC Energy Atlas](https://energy-industry-geolab.jrc.ec.europa.eu/energy-atlas/).
- `retrieve_desnz_electricity_consumption`: subnational electricity consumption for Great Britain from the [Department for Energy Security and Net Zero](https://www.gov.uk/government/statistics/regional-and-local-authority-electricity-consumption-statistics).
- `retrieve_ons_lad`: shapefiles of local authorities in the United Kingdom from the [Office for National Statistics](https://geoportal.statistics.gov.uk/datasets/ons::local-authority-districts-may-2024-boundaries-uk-bsc-2/about).
