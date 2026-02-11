<!--
SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
SPDX-License-Identifier: CC-BY-4.0
-->

![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/open-energy-transition/open-tyndp?include_prereleases)
![Commits behind PyPSA-Eur](https://img.shields.io/github/commits-difference/open-energy-transition/open-tyndp?base=master&head=PyPSA:pypsa-eur:master&label=commits%20behind%20PyPSA-Eur)
[![Test workflows](https://github.com/open-energy-transition/open-tyndp/actions/workflows/test.yaml/badge.svg)](https://github.com/open-energy-transition/open-tyndp/actions/workflows/test.yaml)
[![Documentation](https://readthedocs.org/projects/open-tyndp/badge/?version=latest)](https://open-tyndp.readthedocs.io/en/latest/?badge=latest)
![Size](https://img.shields.io/github/repo-size/open-energy-transition/open-tyndp)
[![Zenodo Open-TYNDP](https://zenodo.org/badge/DOI/10.5281/zenodo.18494363.svg)](https://doi.org/10.5281/zenodo.18494363)
[![Zenodo PyPSA-Eur](https://zenodo.org/badge/DOI/10.5281/zenodo.3520874.svg)](https://doi.org/10.5281/zenodo.3520874)
[![Snakemake](https://img.shields.io/badge/snakemake-≥9-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![Zenodo TYNDP data](https://zenodo.org/badge/DOI/10.5281/zenodo.14230568.svg)](https://doi.org/10.5281/zenodo.14230568)
[![Discord](https://img.shields.io/discord/911692131440148490?logo=discord)](https://discord.gg/AnuJBk23FU)
[![REUSE status](https://api.reuse.software/badge/github.com/open-energy-transition/open-tyndp)](https://api.reuse.software/info/github.com/open-energy-transition/open-tyndp)

# Open-TYNDP: Interfacing Open Energy System Planning with ENTSO-E Models and Contributing to TYNDP
<img src="https://raw.githubusercontent.com/open-energy-transition/oet-website/main/assets/img/oet-logo-red-n-subtitle.png" alt="Open Energy Transition Logo" width="260" height="100" align="right">

> [!WARNING]
> Open-TYNDP is under active development and is not yet feature-complete. The current [development status](#development-status) and the general [Limitations](https://open-tyndp.readthedocs.io/en/latest/limitations.html) are important to understand before using the model.

This repository introduces the open model dataset of the Open-TYNDP research and innovation project, which is a collaboration between [Open Energy Transition (OET)](https://openenergytransition.org/) and the European Network of Transmission System Operators for Electricity (ENTSO-E). The project’s aim is to explore and consider the adoption of PyPSA in the Ten-Year Network Development Plan (TYNDP) by building a workflow based on [PyPSA-Eur](https://github.com/pypsa/pypsa-eur). It seeks to complement the tools currently used in the TYNDP cycles, especially for Scenario Building (SB) and Cost-Benefit Analysis (CBA). This approach is designed to enhance transparency and lower barriers to stakeholder participation in European energy planning. Beyond Europe, the project aspires to demonstrate the viability of open-source (OS) frameworks in energy planning, encouraging broader global adoption.

To build trust in and ensure reproducibility with the new open-source toolchain, the project first focuses on replicating key figures from the 2024 TYNDP cycle, before aligning with the current 2026 TYNDP cycle. This process involves developing new features within the open-source domain to address existing gaps, integrating tools for data interoperability and dynamic visualizations, and publishing best practices to encourage the adoption of open energy models. Additionally, the project emphasizes stakeholder consultations and [interactive workshops](https://open-energy-transition.github.io/open-tyndp-workshops/intro.html) alongside the development of the PyPSA tool, further promoting collaboration and transparency throughout the process.

This repository is a soft-fork of [OET/PyPSA-Eur](https://github.com/open-energy-transition/pypsa-eur) and contains the entire project `Open-TYNDP` supported by OET, including code and documentation. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

This repository is maintained using [OET's soft-fork strategy](https://open-energy-transition.github.io/handbook/docs/Engineering/SoftForkStrategy). OET's primary aim is to contribute as much as possible to the open-source upstream repositories. For long-term changes that cannot be directly merged upstream, the strategy organizes and maintains OET forks, ensuring they remain up-to-date and compatible with upstream on a regular basis, while also supporting future contributions back to the OS repositories.

# Development status

**Warning**: Open-TYNDP is under active development and is not yet feature-complete. The current development status and general [Limitations](https://open-tyndp.readthedocs.io/en/latest/limitations.html) are important to understand before using the model. The model includes partial data from the TYNDP 2024 cycle, and its validation is ongoing. The github repository [issues](https://github.com/open-energy-transition/open-tyndp/issues) collect known topics we are working on (please feel free to help or make suggestions). The fact that this project relies on a soft-fork strategy implies that [upstream issues](https://github.com/PyPSA/PyPSA-Eur/issues) need to be addressed in the PyPSA-Eur repository. The [documentation](https://open-tyndp.readthedocs.io/) also remains work in progress.

## ✅ Already implemented features

The back-casting of the 2024 TYNDP cycle involves developing new features based on the published [modelling methodology report](https://2024.entsos-tyndp-scenarios.eu/wp-content/uploads/2025/01/TYNDP_2024_Scenarios_Methodology_Report_Final_Version_250128.pdf). Major and already implemented features are summarized below. Please, refer to the [Release Notes](https://open-tyndp.readthedocs.io/en/latest/release_notes.html) for a more comprehensive list of features and to the relevant [pull requests](https://github.com/open-energy-transition/open-tyndp/pulls?q=is%3Apr+label%3A%22major+feature%22) for extensive documentation of the implementations.

- Introduced a benchmarking framework that assesses Open-TYNDP quality against published TYNDP 2024 data (see PR [#73](https://github.com/open-energy-transition/open-tyndp/pull/73)) and Visualisation Platform data (see PRs [#117](https://github.com/open-energy-transition/open-tyndp/pull/117), [#281](https://github.com/open-energy-transition/open-tyndp/pull/281)).
- Introduced the electricity and hydrogen reference grids from TYNDP 2024 grid data (see PRs [#18](https://github.com/open-energy-transition/open-tyndp/pull/18), [#17](https://github.com/open-energy-transition/open-tyndp/pull/17), [#340](https://github.com/open-energy-transition/open-tyndp/pull/340)), including hydrogen import potentials and corridors from outside of the modelled countries (see PR [#36](https://github.com/open-energy-transition/open-tyndp/pull/36)). 
- Introduced TYNDP offshore wind hubs topology with both electric and hydrogen infrastructure, offshore electrolysers, and detailed wind farm characteristics (see PR [#54](https://github.com/open-energy-transition/open-tyndp/pull/54)).
- Added processing, preparation and attaching of PECD v3.1 renewable profiles and PEMMDB v2.4 hydro inflows to their corresponding components (see PRs [#53](https://github.com/open-energy-transition/open-tyndp/pull/53), [#71](https://github.com/open-energy-transition/open-tyndp/pull/71), [#77](https://github.com/open-energy-transition/open-tyndp/pull/77), [#338](https://github.com/open-energy-transition/open-tyndp/pull/338)).
- Added processing and preparation of PEMMDB v2.4 capacity, must-run, and availability data, along with expansion trajectories for conventional and renewable power generation, electrolysers, batteries, and DSR (see PR [#97](https://github.com/open-energy-transition/open-tyndp/pull/97)).
- Added the TYNDP electricity, hydrogen, methane demands as exogenous demands (see PRs [#14](https://github.com/open-energy-transition/open-tyndp/pull/14) for electricity; [#169](https://github.com/open-energy-transition/open-tyndp/pull/169) and [#230](https://github.com/open-energy-transition/open-tyndp/pull/230) for hydrogen ; and [#208](https://github.com/open-energy-transition/open-tyndp/pull/208) and [#220](https://github.com/open-energy-transition/open-tyndp/pull/220) for methane).
- Added option to use the TYNDP H2 topology including the TYNDP H2 reference grid, H2 Z1 and Z2 setup, production, reconversion and storage technologies (see PR [#17](https://github.com/open-energy-transition/open-tyndp/pull/17)).
- Added TYNDP solar, onwind, hydro and conventional thermal generation using PEMMDB capacities, must-runs and availabilities (see PRs [#115](https://github.com/open-energy-transition/open-tyndp/pull/115) and [#139](https://github.com/open-energy-transition/open-tyndp/pull/139) for solar and onwind ; [#338](https://github.com/open-energy-transition/open-tyndp/pull/338) for hydro ; [#195](https://github.com/open-energy-transition/open-tyndp/pull/195) for conventionals) and electrolyzer PEMMDB capacities (only for NT scenario) (see PR [#248](https://github.com/open-energy-transition/open-tyndp/pull/248)).
- Added CO2 emission prices per planning horizon (see PR [#198](https://github.com/open-energy-transition/open-tyndp/pull/198)).
- Introduced a workflow structure for performing Cost-Benefit Analysis (CBA) using both TOOT (Take One Out at a Time) and PINT (Put In at a Time) methodologies for TYNDP transmission and storage projects on top of the SB results (see PR [#149](https://github.com/open-energy-transition/open-tyndp/pull/149)).
- Added the PINT and TOOT reference and project network preparation for CBA (see PRs [#199](https://github.com/open-energy-transition/open-tyndp/pull/199), [#211](https://github.com/open-energy-transition/open-tyndp/pull/211), [#353](https://github.com/open-energy-transition/open-tyndp/pull/353)).
- Added weekly rolling horizon optimization for CBA networks (see PR [#217](https://github.com/open-energy-transition/open-tyndp/pull/217), [#385](https://github.com/open-energy-transition/open-tyndp/pull/385)).
- Introduced calculation of B1, B2, B3 and B4 indicator in the CBA (see PRs [#186](https://github.com/open-energy-transition/open-tyndp/pull/186), [#348](https://github.com/open-energy-transition/open-tyndp/pull/348), [#350](https://github.com/open-energy-transition/open-tyndp/pull/350), [#398](https://github.com/open-energy-transition/open-tyndp/pull/398)).
- Added automated Windows installer for easy setup on Windows systems (see PR [#333](https://github.com/open-energy-transition/open-tyndp/pull/333)).
- Automatically launched PyPSA-Explorer to investigate run results (see PR [#443](https://github.com/open-energy-transition/open-tyndp/pull/443)).



|       Feature        | TYNDP 2024 topology | Open-TYNDP topology |
|:--------------------:|:-------------------:|:-------------------:|
| **Electricity Grid** | <img src="doc/img/tyndp/electricity-grid-report.png" height="300px" alt="TYNDP 2024 electricity topology"> | <img src="doc/img/tyndp/electricity-grid.png" height="300px" alt="Open-TYNDP electricity topology"> |
|  **Hydrogen Grid**   | <img src="doc/img/tyndp/h2-grid-report.png" height="300px" alt="TYNDP 2024 hydrogen topology"> | <img src="doc/img/tyndp/h2-grid.png" height="300px" alt="Open-TYNDP hydrogen topology"> |
|  **Offshore Grid**   | <img src="doc/img/tyndp/offshore-grid-report.png" height="300px" alt="TYNDP 2024 offshore topology"> | <img src="doc/img/tyndp/offshore-grid.png" height="300px" alt="Open-TYNDP offshore topology"> |

## 🔨 Features in development

While multiple TYNDP features are already introduced to the Open-TYNDP model, there are several other features and assumptions that are still in development and currently rely on default implementations and assumptions from PyPSA-Eur.

|                       Milestone                        |                       Feature                       |                                                                                                                                             PR                                                                                                                                             | Status |
|:------------------------------------------------------:|:---------------------------------------------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:------:|
|       **Visualizations and workflow automation**       |                 Automated workflow                  |                                                                                                                                                                                                                                                                                            |   ✅    |
|                                                        |               TYNDP plotting routines               |                                                                                            [#443](https://github.com/open-energy-transition/open-tyndp/pull/443)                                                                                                                                                                                                  |   🔨   |
|           **Automated tests and benchmarks**           |           Automated benchmarking routine            |                                     [#73](https://github.com/open-energy-transition/open-tyndp/pull/73), [#117](https://github.com/open-energy-transition/open-tyndp/pull/117), [#281](https://github.com/open-energy-transition/open-tyndp/pull/281)                                      |   ✅    |
|              **TYNDP modelling features**              |           Perfect foresight optimization            |                                                                                                                                                                                                                                                                                            |   ⌛    |
|                                                        |            Security of Supply (SoS) loop            |                                                                                                                                                                                                                                                                                            |   ⌛    |
| **Existing infrastructure and associated parameters**  |             Electricity reference grid              |                                                                         [#18](https://github.com/open-energy-transition/open-tyndp/pull/18), [#340](https://github.com/open-energy-transition/open-tyndp/pull/340)                                                                         |   ✅    |
|                                                        |               Hydrogen reference grid               |                                                                          [#17](https://github.com/open-energy-transition/open-tyndp/pull/17), [#36](https://github.com/open-energy-transition/open-tyndp/pull/36)                                                                          |   ✅    |
|                                                        |                    Offshore grid                    |                                                                                                            [#54](https://github.com/open-energy-transition/open-tyndp/pull/54)                                                                                                             |   ✅    |
|                                                        |                      PECD data                      |                                                                          [#53](https://github.com/open-energy-transition/open-tyndp/pull/53), [#71](https://github.com/open-energy-transition/open-tyndp/pull/71)                                                                          |   ✅    |
|                                                        |                    Hydro inflows                    |                                                                                                            [#77](https://github.com/open-energy-transition/open-tyndp/pull/77), [#338](https://github.com/open-energy-transition/open-tyndp/pull/338)                                                                                                             |   ✅    |
|                                                        |      PEMMDB capacities & must-runs processing       |                                                                                                            [#97](https://github.com/open-energy-transition/open-tyndp/pull/97)                                                                                                             |   ✅    |
|                                                        |    Investment candidates trajectories processing    |                                                                                                            [#97](https://github.com/open-energy-transition/open-tyndp/pull/97)                                                                                                             |   ✅    |
|                    **TYNDP demand**                    |                     Electricity                     |                                                                                                            [#14](https://github.com/open-energy-transition/open-tyndp/pull/14)                                                                                                             |   ✅    |
|                                                        |                      Hydrogen                       |                                                                        [#169](https://github.com/open-energy-transition/open-tyndp/pull/169), [#230](https://github.com/open-energy-transition/open-tyndp/pull/230)                                                                        |   ✅    |
|                                                        |                       Methane                       |                                                                        [#208](https://github.com/open-energy-transition/open-tyndp/pull/208), [#220](https://github.com/open-energy-transition/open-tyndp/pull/220)                                                                        |   🔨    |
|                                                        |                   Synthetic fuels                   |                                                                                                                                                                                                                                                                                            |   ⌛    |
|                                                        |                  District heating                   |                                                                                                                                                                                                                                                                                            |   ⌛    |
|                                                        |                   Energy imports                    |                                                                                                                                                                                                                                                                                            |   ⌛    |
|          **TYNDP technologies and carriers**           |                  TYNDP generation technologies (incl. SRES and DRES)                  | [#115](https://github.com/open-energy-transition/open-tyndp/pull/115), [#139](https://github.com/open-energy-transition/open-tyndp/pull/139), [#195](https://github.com/open-energy-transition/open-tyndp/pull/195), [#248](https://github.com/open-energy-transition/open-tyndp/pull/248), [#338](https://github.com/open-energy-transition/open-tyndp/pull/338), [#447](https://github.com/open-energy-transition/open-tyndp/pull/447) |   🔨   |
|                                                        |                 Prosumer modelling                  |                                                                                                                                                                                                                                                                                            |   ⌛    |
|                                                        |                    EV modelling                     |                                                                                                                                                                                                                                                                                            |   ⌛    |
|                                                        |               Synthetic fuel carriers               |                                                                                                                                                                                                                                                                                            |   ⌛    |
|                                                        |                  Hybrid heat pumps                  |                                                                                                                                                                                                                                                                                            |   ⌛    |
|                                                        |     Industrial electricity and hydrogen demands     |                                                                                                                                                                                                                                                                                            |   ⌛    |
|                                                        |                   Hydrogen zones                    |                                                                                                            [#17](https://github.com/open-energy-transition/open-tyndp/pull/17)                                                                                                             |   ✅    |
|            **CBA Assessment Framework**                |                TOOT/PINT methodology                          | [#149](https://github.com/open-energy-transition/open-tyndp/pull/149), [#199](https://github.com/open-energy-transition/open-tyndp/pull/199), [#211](https://github.com/open-energy-transition/open-tyndp/pull/211), [#217](https://github.com/open-energy-transition/open-tyndp/pull/217), [#353](https://github.com/open-energy-transition/open-tyndp/pull/353), [#385](https://github.com/open-energy-transition/open-tyndp/pull/385), [#426](https://github.com/open-energy-transition/open-tyndp/pull/426) |   🔨    |
|                                                        |             Climate years (weighted average)                  |                                                                                                                                                                                                                                                                                            |   🔨    |
|                                                        |                   CBA Benchmarking                            |                                                                                           [#405](https://github.com/open-energy-transition/open-tyndp/pull/405)                                                                                                                            |   ✅    |
|            **CBA Benefit Indicators**                  |       B1: Socioeconomic Welfare (SEW)                         |                                                                                                           [#186](https://github.com/open-energy-transition/open-tyndp/pull/186)                                                                                                            |   ✅    |
|                                                        |       B2: CO₂ Variation with societal costs                   |                                                                                                                                                                                  [#348](https://github.com/open-energy-transition/open-tyndp/pull/348)                                                                                                                                                                                   |   ✅    |
|                                                        |       B3: DRES Integration (curtailment reduction)            |                                                                                                                                                                                  [#350](https://github.com/open-energy-transition/open-tyndp/pull/350)                                                                                                                                                                                   |   ✅    |
|                                                        |       B4: Non-CO₂ Emissions (NOx, SOx, PM2.5/10, NMVOC, NH₃)  |                                                                                                                                                                                  [#398](https://github.com/open-energy-transition/open-tyndp/pull/398)                                                                                                                                                                                   |   ✅    |
|                                                        |       B5: Grid Losses, B6: Adequacy, B7: Flexibility, B8: Stability, B9: Reserves reduction |                                                                                                                                                                                                                                                                                            |        |

**See also**: See also the [GitHub repository issues](https://github.com/open-energy-transition/open-tyndp/issues) for a comprehensive list of currently open issues.

# Repository structure

* `benchmarks`: will store `snakemake` benchmarks (does not exist initially)
* `config`: configurations used in the study
* `cutouts`: will store raw weather data cutouts from `atlite` (does not exist initially)
* `data`: includes input data that is not produced by any `snakemake` rule
* `doc`: includes all files necessary to build the `readthedocs` documentation of PyPSA-Eur
* `envs`: includes backup `conda` environments if `pixi` installation does not work.
* `logs`: will store log files (does not exist initially)
* `notebooks`: includes all the `notebooks` used for ad-hoc analysis
* `report`: contains all files necessary to build the report; plots and result files are generated automatically
* `rules`: includes all the `snakemake`rules loaded in the `Snakefile`
* `resources`: will store intermediate results of the workflow which can be picked up again by subsequent rules (does not exist initially)
* `results`: will store the solved PyPSA network data, summary files and plots (does not exist initially)
* `scripts`: includes all the Python scripts executed by the `snakemake` rules to build the model

# Installation and Usage

## 1. Installation

### Option A: Windows Installer (Recommended for Windows)

For Windows users, download the latest installer (e.g., `open-tyndp-0.4.0-pixi-Windows-x86_64.exe`) from the [releases page](https://github.com/open-energy-transition/open-tyndp/releases) and run it. The installer automatically sets up everything you need, including pixi, the repository, and the conda environment. See `utils/windows-installer/README.md` for details.

### Option B: Manual Installation (All Platforms)

Clone the repository:

```sh
git clone https://github.com/open-energy-transition/open-tyndp
```
PyPSA-Eur, and consequently Open-TYNDP, relies on a set of other Python packages to function.
We manage these using [pixi](https://pixi.sh/latest/>).
Once pixi is installed, you can activate the project environment for your operating system and have access to all the PyPSA-Eur dependencies from the command line:

```sh
pixi shell -e open-tyndp
```

>[!NOTE]
>`pixi` will create a distinct environment in every project directory, even if you have identical copies of a project cloned locally.
>As there is a common system-level package cache, `pixi` efficiently conserves disk space in such cases.

>[!TIP]
>If `pixi` isn't working, you can install from one of the fallback `conda` environment files found in `envs`.
>For more details see [the PyPSA-Eur installation guide](https://pypsa-eur.readthedocs.io/en/latest/installation.html).

## 2. Run the analysis

```sh
pixi run tyndp
```

This will run all analysis steps to reproduce results and build the report.

To list all the rules that need to be executed (dry run), run:

```sh
pixi run tyndp -n
```

Note: The workflow automatically generates dependency graphs after successful completion (via `rulegraphs` and `filegraphs` rules). The generated graphs are saved to the `resources/` directory.

# Contributing and Support
We strongly welcome anyone interested in contributing to this project. If you have any ideas, suggestions or encounter problems, feel invited to file issues or make pull requests on GitHub.
-   To **discuss** with other PyPSA users, organise projects, share news, and get in touch with the community you can use the [Discord server](https://discord.gg/AnuJBk23FU). Open-TYNDP has its own dedicated channel [pypsa-open-tyndp](https://discord.com/channels/911692131440148490/1414977512089321564) for project-specific discussions.
-   For **bugs and feature requests**, please use the appropriate GitHub issues page. Issues specific to Open-TYNDP belong on the [Open-TYNDP Issues page](https://github.com/open-energy-transition/open-tyndp/issues), while PyPSA-Eur issues should be submitted to the [PyPSA-Eur Github Issues page](https://github.com/PyPSA/pypsa-eur/issues).

# Contact
For any questions about Open-TYNDP or other queries, reach out via the [pypsa-open-tyndp](https://discord.com/channels/911692131440148490/1414977512089321564) channel or <a href="mailto:tyndp@openenergytransition.org">tyndp@openenergytransition.org</a>.

Sign Up for the [project newsletter](https://openenergytransitionnewsletter.eo.page/tyndp-oet) for updates!

# Citation

If you want to cite a specific Open-TYNDP version, since v0.5, each release is archived on Zenodo with a release-specific DOI:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18494363.svg)](https://doi.org/10.5281/zenodo.18494363)

If you use Open-TYNDP in your research, please cite it as shown on Zenodo or using the "Cite this repository" button in the sidebar of [Open-TYNDP](https://github.com/open-energy-transition/open-tyndp).

This work builds upon [PyPSA-Eur](https://github.com/pypsa/pypsa-eur) and follows the methodology described in [ENTSO-E's and entsog's TYNDP 2024 Scenarios Methodology Report](https://2024.entsos-tyndp-scenarios.eu/scenarios-methodology-report/).

# Licence

Open-TYNDP is a soft-fork of PyPSA-Eur, relying on a similar licensing strategy. As with PyPSA-Eur, the code in Open-TYNDP is released as free software under the
[MIT License](https://opensource.org/licenses/MIT), see [`Licenses`](https://open-tyndp.readthedocs.io/en/latest/licenses.html) for attribution and licensing strategy details.
Additionally, different licenses and terms of use may apply to the various
input data, see [`Data Sources`](https://open-tyndp.readthedocs.io/en/latest/data_sources.html).
