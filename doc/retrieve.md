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

The cache is filled in a separate step, once per workflow:

```console
$ pixi run collect-data          # everything Scenario Building needs
$ pixi run collect-data-cba      # everything the CBA needs
$ pixi run tyndp-sb              # reads the cache, no network required
```

Each task dry-runs its own target to see which `retrieve_*` rules that workflow's graph actually
contains and downloads those, with the cache and retrieval forced on so it fills whatever your
config says. It is possible to dry-run this task with `pixi run collect-data -n` to list what
would be downloaded without fetching anything. 
The source selection still applies: with `data_config: tyndp` the cache fills from the
Google Cloud Storage mirror, into the same paths.

!!! warning "The CBA scope is currently only partly knowable in advance"
    `clean_projects` is a Snakemake checkpoint, so the CBA graph gains dependency on SB jobs only 
    once it has run and `pixi run collect-data-cba` sees just the datasets needed up to that point. 
    Run `pixi run collect-data` first and re-run `collect-data-cba` after.

To reach an offline machine, fill the cache where the network is available and copy the
directory across, e.g. using rsync:

```console
$ rsync -a data/local-cache/ offline-machine:~/open-tyndp/data/local-cache/
```

### Reading the cache without a network {#local_cache_offline}

A cached run contacts no remote URL. The retrieve rules stay defined but take the cache
manifest `{local_cache.directory}/.collected` as their input instead of a URL, because Snakemake
otherwise probes remote inputs over the network while the graph is built; filling the cache also drops
Snakemake's provenance records for the cached files, so that switch does not mark every retrieve
job out of date. The manifest and the provenance reset are written only when a collect task
finishes, so an interrupted `collect-data` leaves stale records behind and the next run will
want to re-fetch.

Should the cache be missing a file on offline execution, the workflow stops at parse time and names the absent
files, also on dry runs. It can only report what the cache was told to hold, though: a cache collected for
Scenario Building and then used for a CBA run passes the check and surfaces the gap later as a connection error,
which `pixi run collect-data-cba` fixes.

Data already retrieved into `data/{dataset}/{source}/{version}` is not reused, so the first
`collect-data` downloads it again unless you move those folders into the cache layout by hand.

## Retrieve rules

Below some specific `retrieve_<dataset>` rules are documented.
For more information on the datasets retrieved, see the [data sources](data_sources.md) and *Data inventory* section there in the documentation.

### Rule `retrieve_bidding_zones`

<!-- ::: retrieve_bidding_zones (module not found) -->

### Rule `retrieve_cutout`

See [cutouts](configuration.md#atlite_cf).


### Rule `retrieve_electricity_demand_energy_atlas`

This rule downloads 1km by 1km raster of estimated annual electricity demand from the [JRC Energy Atlas](https://energy-industry-geolab.jrc.ec.europa.eu/energy-atlas/).

### Rule `retrieve_desnz_electricity_consumption`

This rule downloads subnational electricity consumption data for Great Britain from the [Department for Energy Security and Net Zero](https://www.gov.uk/government/statistics/regional-and-local-authority-electricity-consumption-statistics).

### Rule `retrieve_ons_lad`

This rule downloads shapefiles of local authorities in the United Kingdom from the [Office for National Statistics](https://geoportal.statistics.gov.uk/datasets/ons::local-authority-districts-may-2024-boundaries-uk-bsc-2/about).

### Rule `retrieve_electricity_demand_opsd`

This rule downloads hourly electric load data for each country from the [OPSD platform](https://data.open-power-system-data.org/time_series/2019-06-05/time_series_60min_singleindex.csv).

**Relevant Settings**

None.

**Outputs**

- `data/electricity_demand_opsd_raw.csv`

### Rule `retrieve_electricity_demand_entsoe`

This rule downloads hourly electric load data for each country from the [ENTSOE Transparency Platform](https://transparency.entsoe.eu).

**Relevant Settings**

None.

**Outputs**

- `data/electricity_demand_entsoe_raw.csv`

### Rule `retrieve_electricity_demand_neso`

This rule downloads hourly electric load data for the United Kingdom from the [NESO Data Portal](https://www.neso.energy/data-portal/historic-demand-data).

**Relevant Settings**

None.

**Outputs**

- `data/electricity_demand_neso_raw.csv`


### Rule `retrieve_cost_data`

This rule downloads techno-economic assumptions from the [technology-data repository](https://github.com/pypsa/technology-data).

**Relevant Settings**

```yaml
costs:
    year:
```

!!! seealso
    Documentation of the configuration file `config/config.yaml` at
    [costs_cf](configuration.md#costs_cf)

**Outputs**

- `data/costs/primary/{version}/costs_{year}.csv`
