..
  SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>

  SPDX-License-Identifier: CC-BY-4.0

.. _installation:

##########################################
Installation
##########################################

The subsequently described installation steps are demonstrated as shell commands, where the path before the ``%`` sign denotes the
directory in which the commands following the ``%`` should be entered.

Clone the Repository
====================

First of all, clone the `Open-TYNDP repository <https://github.com/open-energy-transition/open-tyndp>`__ using the version control system ``git`` in the command line.

.. code:: console

    $ git clone https://github.com/open-energy-transition/open-tyndp.git


.. _deps:

Install Python Dependencies
===============================

PyPSA-Eur, and consequently Open-TYNDP, relies on a set of other Python packages to function. We recommend
using the package manager `conda <https://docs.anaconda.com/miniconda/>`__ or
`mamba <https://mamba.readthedocs.io/en/latest/>`__ to install them and manage
your environments.

The package requirements are curated in the ``envs/environment.yaml`` file.
There are also regularly updated locked environment files for each platform generated with conda-lock to
ensure reproducibility. Choose the correct file for your platform:

* For Intel/AMD processors:

  - Linux: ``envs/linux-64.lock.yaml``

  - macOS: ``envs/osx-64.lock.yaml``

  - Windows: ``envs/win-64.lock.yaml``

* For ARM processors:

  - macOS (Apple Silicon): ``envs/osx-arm64.lock.yaml``

  - Linux (ARM): Currently not supported via lock files; requires building certain packages, such as ``PySCIPOpt``, from source

We recommend using these locked files for a stable environment.

.. code:: console

    $ conda update conda

    $ conda env create -n open-tyndp -f envs/linux-64.lock.yaml # select the appropriate file for your platform

    $ conda activate open-tyndp

Install a Solver
================

PyPSA passes the PyPSA-Eur network model to an external solver for performing the optimisation.
PyPSA is known to work with the free software

- `HiGHS <https://highs.dev/>`__
- `Cbc <https://projects.coin-or.org/Cbc#DownloadandInstall>`__
- `GLPK <https://www.gnu.org/software/glpk/>`__ (`WinGLKP <http://winglpk.sourceforge.net/>`__)
- `SCIP <https://scipopt.github.io/PySCIPOpt/docs/html/index.html>`__

and the non-free, commercial software (for some of which free academic licenses are available)

- `Gurobi <https://www.gurobi.com/documentation/quickstart.html>`__
- `CPLEX <https://www.ibm.com/products/ilog-cplex-optimization-studio>`__
- `FICO Xpress Solver <https://www.fico.com/de/products/fico-xpress-solver>`__

For installation instructions of these solvers for your operating system, follow the links above.
Commercial solvers such as Gurobi and CPLEX currently significantly outperform open-source solvers for large-scale problems, and
it might be the case that you can only retrieve solutions by using a commercial solver.
Nevertheless, you can still use open-source solvers for smaller problems.

.. seealso::
    `Instructions how to install a solver in the documentation of PyPSA <https://pypsa.readthedocs.io/en/latest/installation.html#getting-a-solver-for-linear-optimisation>`__

.. note::
    The rules :mod:`cluster_network` solves a mixed-integer quadratic optimisation problem for clustering.
    The open-source solvers HiGHS, Cbc and GlPK cannot handle this. A fallback to SCIP is implemented in this case, which is included in the standard environment specifications.
    For an open-source solver setup install for example HiGHS **and** SCIP in your ``conda`` environment on OSX/Linux.
    To install the default solver Gurobi, run

    .. code:: console

        $ conda activate pypsa-eur
        $ conda install -c gurobi gurobi"=12.0.1"

    Additionally, you need to setup your `Gurobi license <https://www.gurobi.com/solutions/licensing/>`__.
