<!-- SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# TYNDP 2026 data

The TYNDP 2026 data is downloaded and extracted using the `retrieve_tyndp_2026` and `retrieve_tyndp_pecd` rules. Originally published by ENTSO-E and ENTSOG under
Creative Commons Attribution 4.0 International License (CC-BY 4.0).
The original data files can be found under <https://2026.entsos-tyndp-scenarios.eu/download/>.
The `{version}`, `{pecd_version}` and `{pecd_prebuilt_version}` placeholders refer to the dataset versions listed in `data/tyndp_versions.csv`.

## `Demand`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/inputs/Demand.zip`
- **Description:** Contains hourly electricity, hydrogen, synthetic fuel, behind-the-meter, EV fleet and thermal (hydrogen and methane) demand profiles by node and weather scenario, for 2030, 2035, 2040 and 2050.

## `Electric-Vehicle`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/inputs/Electric-Vehicle.zip`
- **Description:** Contains the passenger EV fleet inputs: EV stock and demand, market zone, flexible fleet and V2G shares, charging station and battery parameters, and availability, driving and minimum state of charge profiles.

## `Heat`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/inputs/Heat.zip`
- **Description:** Contains coefficient of performance (COP) profiles per planning horizon and the hybrid heat pump parameters (heat pumps and methane or hydrogen boilers).

## `Hydro-Inflow`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/inputs/Hydro-Inflow.zip`
- **Description:** Contains PECD v4.2 hydro inflow time series and hydro constraints per weather scenario, for 2030, 2035, 2040 and 2050.

## `Hydrogen`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/inputs/Hydrogen.zip`
- **Description:** Contains hydrogen storage parameters, steam methane reforming (SMR and SMR+CCS) assets, hydrogen import generators and profiles, and synthetic fuel lines and generator properties.

## `Line-data`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/inputs/Line-data.zip`
- **Description:** Contains the electricity and hydrogen reference grids and the wheeling charges.

## `Nodes`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/inputs/Nodes.zip`
- **Description:** Contains a list of electricity and hydrogen nodes.

## `PEMMDB_2.X`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/inputs/PEMMDB_2.X.zip`
- **Description:** Contains the PEMMDB 2.X installed generation and storage capacities, must-run constraints and technical parameters per node, for 2030, 2035, 2040 and 2050.

## `PECD/PECD_{pecd_version}`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/PECD/PECD_{pecd_version}.zip`
- **Description:** Contains hourly PECD v4.2 capacity factor time series for Solar PV rooftop,
  Solar PV utility, Onshore Wind, Offshore Wind and Solar CSP.

## `PECD/PECD_{pecd_prebuilt_version}`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/PECD/PECD_{pecd_prebuilt_version}.zip`
- **Description:** Contains pre-built hourly PECD v4.2 capacity factor time series for weather scenarios WS003, WS021, WS029, WS032, WS037, WS059, WS065, WS071, WS077, WS091, WS092 and WS106, covering Solar PV rooftop, Solar PV utility, Onshore Wind, Offshore Wind and Solar CSP.

## `Supply-Tool`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/outputs/Supply-Tool.zip`
- **Description:** Contains the TYNDP 2026 Supply Tool, with one workbook per scenario (NT+, HEV and LEV).

## `NT+_TimeSeriesDashboard`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/outputs/NT+_TimeSeriesDashboard.zip`
- **Description:** Contains the hourly market model outputs of the NT+ scenario, one file per country, for each planning horizon and weather scenario. Used for benchmarking and to build the nuclear availability profiles: the hourly nuclear generation of the `E-Market` sheet is divided by the installed capacity of the `Installed Capacity` sheet, and the weekly maximum is retained to account for maintenance schedules and forced outages.

## `TYNDP_2026_Scenario_Report_Data_Figures.xlsx`

- **Source:** TYNDP 2026 Scenarios, <https://2026.entsos-tyndp-scenarios.eu/download/>
- **Link:** `https://storage.googleapis.com/open-tyndp-data-store/2026/{version}/additional/TYNDP_2026_Scenario_Report_Data_Figures.xlsx`
- **Description:** Contains the data used to produce the figures in the TYNDP 2026 Scenarios Report. Used for benchmarking.
