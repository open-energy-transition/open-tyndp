<!-- SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Custom CBA projects

Templates for evaluating user-defined projects with the CBA workflow. Each file ships with its
header only; add one row per project, bus or generator to include it in a run.

| File | Purpose |
| --- | --- |
| `transmission_projects.csv` | Modify an existing transmission project, or add a new one |
| `buses.csv` | Buses added for a custom project, for components that do not connect to an existing bus |
| `generators_static.csv` | Static attributes of custom generators associated with a transmission or storage project |
| `generators_dynamic.csv` | Time-varying attributes of the custom generators added in `generators_static.csv` (two-row header, snapshots as index) |

Custom buses are added before the custom generators, so a generator can be attached to a bus
defined in `buses.csv` within the same run.

Columns, validation rules and defaults are documented under
[Evaluation of custom projects](../../../doc/cba.md#evaluation-of-custom-projects); which projects
are evaluated is configured by `cba.projects`, see
[Selecting custom projects](../../../doc/cba.md#selecting-custom-projects).

The files are read and cleaned by [`scripts/cba/clean_projects.py`](../../../scripts/cba/clean_projects.py),
whose module docstring describes the resulting tables in `resources/cba/`. One can optionally provide the same files as `.xlsx` format for each input file; if detected this will take precedence over given csv files.

Filled-in versions of each file are available in [`examples/`](examples), covering a modified
project, a new project, a new bus and a pair of identical projects used to assess two generator
variants separately. They cover the temporal resolution of the test config `config/test/config.tyndp.yaml` and are not read by the workflow.
