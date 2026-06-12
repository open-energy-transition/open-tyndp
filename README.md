<!--
SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
SPDX-License-Identifier: CC-BY-4.0
-->

![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/open-energy-transition/open-tyndp?include_prereleases)
![Commits behind PyPSA-Eur](https://img.shields.io/github/commits-difference/open-energy-transition/open-tyndp?base=master&head=PyPSA:pypsa-eur:master&label=commits%20behind%20PyPSA-Eur)
[![Test workflows](https://github.com/open-energy-transition/open-tyndp/actions/workflows/test.yaml/badge.svg)](https://github.com/open-energy-transition/open-tyndp/actions/workflows/test.yaml)
[![Documentation](https://readthedocs.org/projects/open-tyndp/badge/?version=latest)](https://open-tyndp.readthedocs.io/en/latest/?badge=latest)
[![Website](https://img.shields.io/badge/website-open--tyndp-4a6fa5?logo=internet-explorer)](https://open-tyndp.openenergytransition.org/)
[![Snakemake](https://img.shields.io/badge/snakemake-≥9-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![Discord](https://img.shields.io/discord/911692131440148490?logo=discord)](https://discord.gg/AnuJBk23FU)
[![REUSE status](https://api.reuse.software/badge/github.com/open-energy-transition/open-tyndp)](https://api.reuse.software/info/github.com/open-energy-transition/open-tyndp)

**Citations & Data:**

[![Cite Open-TYNDP](https://img.shields.io/badge/cite-Open--TYNDP-blue)](https://doi.org/10.5281/zenodo.19372053)
[![Cite PyPSA-Eur](https://img.shields.io/badge/cite-PyPSA--Eur-blue)](https://doi.org/10.5281/zenodo.3520874)
[![Cite PyPSA-Eur-Sec](https://img.shields.io/badge/cite-PyPSA--Eur--Sec-blue)](https://doi.org/10.5281/zenodo.3938042)
[![TYNDP input data](https://img.shields.io/badge/data-TYNDP%202024%20inputs-blue)](https://doi.org/10.5281/zenodo.14230568)

# Open-TYNDP: Interfacing Open Energy System Planning with ENTSO-E Models and Contributing to TYNDP
<img src="https://raw.githubusercontent.com/open-energy-transition/oet-website/main/assets/img/oet-logo-red-n-subtitle.png" alt="Open Energy Transition Logo" width="260" height="100" align="right">

> [!WARNING]
> Open-TYNDP is under active development and is not yet feature-complete. The current [development status](#development-status) and the general [Limitations](https://open-tyndp.readthedocs.io/en/latest/limitations.html) are important to understand before using the model.

This repository introduces the open model dataset of the Open-TYNDP research and innovation project, which is a collaboration between [Open Energy Transition (OET)](https://openenergytransition.org/) and the European Network of Transmission System Operators for Electricity (ENTSO-E). The project’s aim is to explore and consider the adoption of PyPSA in the Ten-Year Network Development Plan (TYNDP) by building a workflow based on [PyPSA-Eur](https://github.com/pypsa/pypsa-eur). It seeks to complement the tools currently used in the TYNDP cycles, especially for Scenario Building (SB) and Cost-Benefit Analysis (CBA). This approach is designed to enhance transparency and lower barriers to stakeholder participation in European energy planning. Beyond Europe, the project aspires to demonstrate the viability of open-source (OS) frameworks in energy planning, encouraging broader global adoption.

To build trust in and ensure reproducibility with the new open-source toolchain, the project first focuses on replicating key figures from the 2024 TYNDP cycle, before aligning with the current 2026 TYNDP cycle. This process involves developing new features within the open-source domain to address existing gaps, integrating tools for data interoperability and dynamic visualizations, and publishing best practices to encourage the adoption of open energy models. Additionally, the project emphasizes stakeholder consultations and [interactive workshops](https://open-energy-transition.github.io/open-tyndp-workshops/intro.html) alongside the development of the PyPSA tool, further promoting collaboration and transparency throughout the process.

This repository is a soft-fork of [OET/PyPSA-Eur](https://github.com/open-energy-transition/pypsa-eur) and contains the entire project `Open-TYNDP` supported by OET, including code and documentation. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

This repository is maintained using [OET's soft-fork strategy](https://open-energy-transition.github.io/handbook/docs/Engineering/SoftForkStrategy). OET's primary aim is to contribute as much as possible to the open-source upstream repositories. For long-term changes that cannot be directly merged upstream, the strategy organizes and maintains OET forks, ensuring they remain up-to-date and compatible with upstream on a regular basis, while also supporting future contributions back to the OS repositories.

# Development status

**Warning**: Open-TYNDP is under active development and is not yet feature-complete. The current [development status](https://open-tyndp.readthedocs.io/en/latest/#development-status) and general [Limitations](https://open-tyndp.readthedocs.io/en/latest/limitations.html) are important to understand before using the model. The model includes partial data from the TYNDP 2024 cycle, and its validation is ongoing. The github repository [issues](https://github.com/open-energy-transition/open-tyndp/issues) collect known topics we are working on (please feel free to help or make suggestions). The fact that this project relies on a soft-fork strategy implies that [upstream issues](https://github.com/PyPSA/PyPSA-Eur/issues) need to be addressed in the PyPSA-Eur repository. The [documentation](https://open-tyndp.readthedocs.io/) also remains work in progress.


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

To run all analysis steps of the Scenario Building to reproduce results and build the report, you can execute:

```sh
pixi run tyndp-sb
```

To list all the rules that need to be executed (dry run), run:

```sh
pixi run tyndp-sb -n
```

Similarly, to run all steps of the Cost-Benefit Analysis, you can execute:

```sh
pixi run tyndp-cba
```

and append `-n` for the corresponding dry run:

```sh
pixi run tyndp-cba -n
```

Note: The workflow automatically generates dependency graphs after successful completion (via `rulegraphs` and `filegraphs` rules). The generated graphs are saved to the `resources/` directory.

# Contributing and Support
We strongly welcome anyone interested in contributing to this project. If you have any ideas, suggestions or encounter problems, feel invited to file issues or make pull requests on GitHub.
-   To **discuss** with other PyPSA users, organise projects, share news, and get in touch with the community you can use the [Discord server](https://discord.gg/AnuJBk23FU). Open-TYNDP has its own dedicated channel [pypsa-open-tyndp](https://discord.com/channels/911692131440148490/1414977512089321564) for project-specific discussions.
-   For **bugs and feature requests**, please [open a new issue](https://github.com/open-energy-transition/open-tyndp/issues/new/choose). We provide templates geared for community members to report bugs or request features, but you're also welcome to use the more detailed maintainer-oriented templates if you're comfortable with them. The Open-TYNDP team reviews new issues shortly after they are submitted. Issues specific to Open-TYNDP belong on the [Open-TYNDP Issues page](https://github.com/open-energy-transition/open-tyndp/issues), while PyPSA-Eur issues should be submitted to the [PyPSA-Eur Github Issues page](https://github.com/PyPSA/pypsa-eur/issues).

# Contact
For any questions about Open-TYNDP or other queries, reach out via the [pypsa-open-tyndp](https://discord.com/channels/911692131440148490/1414977512089321564) channel or <a href="mailto:tyndp@openenergytransition.org">tyndp@openenergytransition.org</a>.

Sign Up for the [project newsletter](https://openenergytransitionnewsletter.eo.page/tyndp-oet) for updates!

# Citation

If you want to cite a specific Open-TYNDP version, since v0.5, each release is archived on Zenodo with a release-specific DOI:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19372053.svg)](https://doi.org/10.5281/zenodo.19372053)

Versions v0.5 and v0.5.1 are archived at [10.5281/zenodo.18494362](https://doi.org/10.5281/zenodo.18494362).

If you use Open-TYNDP in your research, please cite it as shown on Zenodo or using the "Cite this repository" button in the sidebar of [Open-TYNDP](https://github.com/open-energy-transition/open-tyndp).

This work builds upon [PyPSA-Eur](https://github.com/pypsa/pypsa-eur) and follows the methodology described in [ENTSO-E's and ENTSOG's TYNDP 2024 Scenarios Methodology Report](https://2024.entsos-tyndp-scenarios.eu/scenarios-methodology-report/).

# Licence

Open-TYNDP is a soft-fork of PyPSA-Eur, relying on a similar licensing strategy. As with PyPSA-Eur, the code in Open-TYNDP is released as free software under the
[MIT License](https://opensource.org/licenses/MIT), see [`Licenses`](https://open-tyndp.readthedocs.io/en/latest/licenses.html) for attribution and licensing strategy details.
Additionally, different licenses and terms of use may apply to the various
input data, see [`Data Sources`](https://open-tyndp.readthedocs.io/en/latest/data_sources.html).
