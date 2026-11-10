<!-- SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# TYNDP 2026 data

The TYNDP 2026 data is downloaded and extracted using the `retrieve_tyndp_2026` and  `retrieve_tyndp_pecd` rules. Originally published by ENTSO-E and ENTSOG under
Creative Commons Attribution 4.0 International License (CC-BY 4.0).
The original data files can be found under <https://2026.entsos-tyndp-scenarios.eu/download/>.


## `PECD/PECD_{version}`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/PECD/PECD_{version}.zip`

- **Description:** Contains hourly PECD v4.2 capacity factor time series for Solar PV rooftop,
  Solar PV utility, Onshore Wind, Offshore Wind and Solar CSP.

## `PECD/PECD_{pre-built_version}`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/PECD/PECD_{pre-built_version}.zip`
- **Description:** Contains pre-built hourly PECD v4.2 capacity factor time series for weather scenarios WS003, WS021, WS029, WS032, WS037, WS059, WS065, WS071, WS077, WS091, WS092 and WS106, covering Solar PV rooftop, Solar PV utility, Onshore Wind, Offshore Wind and Solar CSP.

## `NT+_TimeSeriesDashboard`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/outputs/NT+_TimeSeriesDashboard.zip`
- **Description:** Contains the hourly market model outputs of the NT+ scenario, one file per country, for each planning horizon and weather scenario. Used for benchmarking and to build the nuclear availability profiles: the hourly nuclear generation of the `E-Market` sheet is divided by the installed capacity of the `Installed Capacity` sheet, and the weekly maximum is retained to account for maintenance schedules and forced outages.

## `cutout/{version}/europe-grid`

- **Source:** Open-TYNDP, generated with [atlite](https://atlite.readthedocs.io)
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/cutout/{version}/europe-grid.nc`
- **Description:** Grid-only cutout without weather data, replacing the ERA5/SARAH-3 cutouts.
  Weather-dependent inputs come from TYNDP data; the cutout only provides the grid for
  `build_population_layouts` and `build_clustered_population_layouts`. The grid is identical
  to the PyPSA-Eur `europe-*-sarah3-era5` cutouts (v1.0). Generated with:

  ```python
  import atlite

  cutout = atlite.Cutout(
      "europe-grid.nc",
      module="era5",
      x=slice(-12.0, 40.8),
      y=slice(33.0, 72.0),
      dx=0.3,
      dy=0.3,
      time="2013-01-01",
  )
  cutout.data.to_netcdf("europe-grid.nc")
  ```
