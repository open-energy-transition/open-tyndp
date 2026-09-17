<!-- SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Custom CBA projects

Templates for evaluating user-defined projects with the CBA workflow. Each file ships with its
header only; add one row per project or generator to include it in a run.

| File | Purpose |
| --- | --- |
| `transmission_projects.csv` | Modify an existing PINT transmission project, or add a new one |
| `generators_static.csv` | Static attributes of custom generators grouped with a transmission or storage project |
| `generators_dynamic.csv` | Time-varying attributes of those generators (two-row header, snapshots as index) |

Columns, validation rules and defaults are documented under
[Evaluation of custom projects](../../../doc/cba.md#evaluation-of-custom-projects); which projects
are evaluated is configured by `cba.projects`, see
[Selecting custom projects](../../../doc/cba.md#selecting-custom-projects).

The files are read and cleaned by [`scripts/cba/clean_projects.py`](../../../scripts/cba/clean_projects.py),
whose module docstring describes the resulting tables in `resources/cba/`.

Filled-in versions of each file are available in [`examples/`](examples), covering a modified
project, a new project and a pair of identical projects used to assess two generator variants
separately. They are not read by the workflow.
