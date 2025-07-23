<!--
SPDX-FileCopyrightText: Open Energy Transition gGmbH and contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
SPDX-License-Identifier: CC-BY-4.0
-->

![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/open-energy-transition/open-tyndp?include_prereleases)
[![Test workflows](https://github.com/open-energy-transition/open-tyndp/actions/workflows/test.yaml/badge.svg)](https://github.com/open-energy-transition/open-tyndp/actions/workflows/test.yaml)
[![Documentation](https://readthedocs.org/projects/open-tyndp/badge/?version=latest)](https://open-tyndp.readthedocs.io/en/latest/?badge=latest)
![Size](https://img.shields.io/github/repo-size/open-energy-transition/open-tyndp)
[![Zenodo PyPSA-Eur](https://zenodo.org/badge/DOI/10.5281/zenodo.3520874.svg)](https://doi.org/10.5281/zenodo.3520874)
[![Zenodo PyPSA-Eur-Sec](https://zenodo.org/badge/DOI/10.5281/zenodo.3938042.svg)](https://doi.org/10.5281/zenodo.3938042)
[![Zenodo TYNDP data](https://zenodo.org/badge/DOI/10.5281/zenodo.14230568.svg)](https://doi.org/10.5281/zenodo.14230568)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.14.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![Discord](https://img.shields.io/discord/911692131440148490?logo=discord)](https://discord.gg/AnuJBk23FU)
[![REUSE status](https://api.reuse.software/badge/github.com/open-energy-transition/open-tyndp)](https://api.reuse.software/info/github.com/open-energy-transition/open-tyndp)

# Open-TYNDP
<img src="https://raw.githubusercontent.com/open-energy-transition/oet-website/main/assets/img/oet-logo-red-n-subtitle.png" alt="Open Energy Transition Logo" width="260" height="100" align="right">

This repository introduces the open model dataset of the Open-TYNDP research and innovation project which is a collaboration between [Open Energy Transition (OET)](https://openenergytransition.org/)<sup>*</sup> and the European Network of Transmission System Operators for Electricity (ENTSO-E). The project’s aim is to explore and consider the adoption of PyPSA in the TYNDP by building a workflow based on PyPSA-Eur. It seeks to complement the tools currently used in the TYNDP cycles, especially for Scenario Building (SB) and Cost-Benefit Analysis (CBA). This approach is designed to enhance transparency and lower barriers to stakeholder participation in European energy planning. Beyond Europe, the project aspires to demonstrate the viability of open-source frameworks in energy planning, encouraging broader global adoption.

To build trust in and ensure reproducibility with the new open-source toolchain, the project first focuses on replicating key figures from the 2024 TYNDP cycle, before aligning with the current 2026 TYNDP cycle. This process involves developing new features within the open-source domain to address existing gaps, integrating tools for data interoperability and dynamic visualizations, and publishing best practices to encourage the adoption of open energy models. Additionally, the project emphasizes stakeholder consultations and [interactive workshops](https://open-energy-transition.github.io/open-tyndp-workshops/intro.html) alongside the development of the PyPSA tool, further promoting collaboration and transparency throughout the process.

This repository is a soft-fork of [OET/PyPSA-Eur](https://github.com/open-energy-transition/pypsa-eur) and contains the entire project `Open-TYNDP` supported by OET, including code and documentation. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

This repository is maintained using [OET's soft-fork strategy](https://open-energy-transition.github.io/handbook/docs/Engineering/SoftForkStrategy). OET's primary aim is to contribute as much as possible to the open source (OS) upstream repositories. For long-term changes that cannot be directly merged upstream, the strategy organizes and maintains OET forks, ensuring they remain up-to-date and compatible with upstream on a regular basis, while also supporting future contributions back to the OS repositories.

# Development status
The back-casting of the 2024 TYNDP cycle involves developing new features based on the published [modelling methodology report](https://2024.entsos-tyndp-scenarios.eu/wp-content/uploads/2025/01/TYNDP_2024_Scenarios_Methodology_Report_Final_Version_250128.pdf). Major and already implemented features are summarized below. Please, refer to the [release notes](https://open-tyndp.readthedocs.io/en/latest/release_notes.html) for a more comprehensive list of features and to the relevant [pull requests](https://github.com/open-energy-transition/open-tyndp/pulls?q=is%3Apr+label%3A%22major+feature%22) for extensive documentation of the implementations.

* Introduced a new electricity base network using TYNDP 2024 electricity reference grid data.
* Added option to use the TYNDP H2 topology including the TYNDP H2 reference grid, H2 Z1 and Z2 setup, production, reconversion and storage technologies.
* Added TYNDP hydrogen import potentials and corridors from outside of the modelled countries. 
* Added the TYNDP electricity demand as an exogenously set demand. 
* Added processing and preparation of TYNDP 2024 PECD v3.1 renewable profiles, replacing default ERA5-based profiles processed with Atlite.
* Introduced TYNDP offshore wind hubs topology with both electric and hydrogen infrastructure, offshore electrolysers, and detailed wind farm characteristics.

|       Feature        | TYNDP 2024 topology | Open-TYNDP topology |
|:--------------------:|:-------------------:|:-------------------:|
| **Electricity Grid** | <img src="doc/img/tyndp/electricity-grid-report.png" height="300px" alt="TYNDP 2024 electricity topology"> | <img src="doc/img/tyndp/electricity-grid.png" height="300px" alt="Open-TYNDP electricity topology"> |
|  **Hydrogen Grid**   | <img src="doc/img/tyndp/h2-grid-report.png" height="300px" alt="TYNDP 2024 hydrogen topology"> | <img src="doc/img/tyndp/h2-grid.png" height="300px" alt="Open-TYNDP hydrogen topology"> |
|  **Offshore Grid**   | <img src="doc/img/tyndp/offshore-grid-report.png" height="300px" alt="TYNDP 2024 offshore topology"> | <img src="doc/img/tyndp/offshore-grid.png" height="300px" alt="Open-TYNDP offshore topology"> |

**WARNING**: Open-TYNDP is under active development and has several
[limitations](https://open-tyndp.readthedocs.io/en/latest/limitations.html) which
you should understand before using the model. The github repository
[issues](https://github.com/open-energy-transition/open-tyndp/issues) collect known topics we are
working on (please feel free to help or make suggestions). The
[documentation](https://open-tyndp.readthedocs.io/) remains somewhat patchy.

# Repository structure

* `benchmarks`: will store `snakemake` benchmarks (does not exist initially)
* `config`: configurations used in the study
* `cutouts`: will store raw weather data cutouts from `atlite` (does not exist initially)
* `data`: includes input data that is not produced by any `snakemake` rule
* `doc`: includes all files necessary to build the `readthedocs` documentation of PyPSA-Eur
* `envs`: includes all the `mamba` environment specifications to run the workflow
* `logs`: will store log files (does not exist initially)
* `notebooks`: includes all the `notebooks` used for ad-hoc analysis
* `report`: contains all files necessary to build the report; plots and result files are generated automatically
* `rules`: includes all the `snakemake`rules loaded in the `Snakefile`
* `resources`: will store intermediate results of the workflow which can be picked up again by subsequent rules (does not exist initially)
* `results`: will store the solved PyPSA network data, summary files and plots (does not exist initially)
* `scripts`: includes all the Python scripts executed by the `snakemake` rules to build the model

# Installation and Usage

## 1. Installation

Clone the repository:

    git clone https://github.com/open-energy-transition/open-tyndp

You need a package manager like [mamba](https://mamba.readthedocs.io/en/latest/) to run the analysis. Users may also prefer to use [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) or [conda](https://docs.conda.io/projects/conda/en/latest/index.html). Using `mamba`, you can create an environment for `<your-os>` from within you can run it:

    mamba env create -n open-tyndp -f envs/<your-os>.lock.yaml

Activate the newly created `open-tydnp` environment:

    mamba activate open-tyndp

## 2. Run the analysis

    make tyndp

This will run all analysis steps to reproduce results and build the report.

To generate a PDF of the dependency graph of all steps `resources/dag_rulegraph.pdf` run:

    snakemake -c1 rulegraph --configfile config/config.tyndp.yaml

<sup>*</sup> Open Energy Transition (g)GmbH, Königsallee 52, 95448 Bayreuth, Germany

----

----

# PyPSA-Eur: A Sector-Coupled Open Optimisation Model of the European Energy System

PyPSA-Eur is an open model dataset of the European energy system at the
transmission network level that covers the full ENTSO-E area. The model is suitable both for operational studies and generation and transmission expansion planning studies.
The continental scope and highly resolved spatial scale enables a proper description of the long-range
smoothing effects for renewable power generation and their varying resource availability.




The model is described in the [documentation](https://pypsa-eur.readthedocs.io)
and in the paper
[PyPSA-Eur: An Open Optimisation Model of the European Transmission
System](https://arxiv.org/abs/1806.01613), 2018,
[arXiv:1806.01613](https://arxiv.org/abs/1806.01613).
The model building routines are defined through a snakemake workflow.
Please see the [documentation](https://pypsa-eur.readthedocs.io/)
for installation instructions and other useful information about the snakemake workflow.
The model is designed to be imported into the open toolbox
[PyPSA](https://github.com/PyPSA/PyPSA).

**WARNING**: PyPSA-Eur is under active development and has several
[limitations](https://pypsa-eur.readthedocs.io/en/latest/limitations.html) which
you should understand before using the model. The github repository
[issues](https://github.com/PyPSA/pypsa-eur/issues) collect known topics we are
working on (please feel free to help or make suggestions). The
[documentation](https://pypsa-eur.readthedocs.io/) remains somewhat patchy. You
can find showcases of the model's capabilities in the Joule paper [The potential
role of a hydrogen network in
Europe](https://doi.org/10.1016/j.joule.2023.06.016), another [paper in Joule
with a description of the industry
sector](https://doi.org/10.1016/j.joule.2022.04.016), or in [a 2021 presentation
at EMP-E](https://nworbmot.org/energy/brown-empe.pdf). We do not recommend to
use the full resolution network model for simulations. At high granularity the
assignment of loads and generators to the nearest network node may not be a
correct assumption, depending on the topology of the underlying distribution
grid, and local grid bottlenecks may cause unrealistic load-shedding or
generator curtailment. We recommend to cluster the network to a couple of
hundred nodes to remove these local inconsistencies. See the discussion in
Section 3.4 "Model validation" of the paper.


![PyPSA-Eur Grid Model](doc/img/elec.png)

The dataset consists of:

- A grid model based on a modified [GridKit](https://github.com/bdw/GridKit)
  extraction of the [ENTSO-E Transmission System
  Map](https://www.entsoe.eu/data/map/). The grid model contains 7072 lines
  (alternating current lines at and above 220kV voltage level and all high
  voltage direct current lines) and 3803 substations.
- The open power plant database
  [powerplantmatching](https://github.com/PyPSA/powerplantmatching).
- Electrical demand time series from the
  [OPSD project](https://open-power-system-data.org/).
- Renewable time series based on ERA5 and SARAH, assembled using the [atlite tool](https://github.com/PyPSA/atlite).
- Geographical potentials for wind and solar generators based on land use (CORINE) and excluding nature reserves (Natura2000) are computed with the [atlite library](https://github.com/PyPSA/atlite).

A sector-coupled extension adds demand
and supply for the following sectors: transport, space and water
heating, biomass, industry and industrial feedstocks, agriculture,
forestry and fishing. This completes the energy system and includes
all greenhouse gas emitters except waste management and land use.

This diagram gives an overview of the sectors and the links between
them:

![sector diagram](doc/img/multisector_figure.png)

Each of these sectors is built up on the transmission network nodes
from [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur):

![network diagram](https://github.com/PyPSA/pypsa-eur/blob/master/doc/img/base.png?raw=true)

For computational reasons the model is usually clustered down
to 50-200 nodes.

Already-built versions of the model can be found in the accompanying [Zenodo
repository](https://doi.org/10.5281/zenodo.3601881).

# Contributing and Support
We strongly welcome anyone interested in contributing to this project. If you have any ideas, suggestions or encounter problems, feel invited to file issues or make pull requests on GitHub.
-   To **discuss** with other PyPSA users, organise projects, share news, and get in touch with the community you can use the [Discord server](https://discord.gg/AnuJBk23FU).
-   For **bugs and feature requests**, please use the [PyPSA-Eur Github Issues page](https://github.com/PyPSA/pypsa-eur/issues).

# Licence

The code in PyPSA-Eur is released as free software under the
[MIT License](https://opensource.org/licenses/MIT), see [`doc/licenses.rst`](doc/licenses.rst).
However, different licenses and terms of use may apply to the various
input data, see [`doc/data_sources.rst`](doc/data_sources.rst).
