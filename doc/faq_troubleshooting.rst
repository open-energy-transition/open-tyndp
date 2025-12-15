.. SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _faq_troubleshooting:

###############################################
FAQ and Troubleshooting
###############################################

This page contains frequently asked questions about the Open-TYNDP project:

General Questions
==================

.. admonition:: What is Open-TYNDP?

   Open-TYNDP is an open-source research and innovation project, which is a collaboration between `Open Energy Transition (OET) <https://openenergytransition.org/>`__ and ENTSO-E. The project aims to explore the option of a complementary open-source tool in the Ten-Year Network Development Plan (TYNDP) by building a workflow based on PyPSA-Eur. It provides one streamlined tool for the Scenario Building (SB) and the Cost-Benefit Analysis (CBA) of the TYNDP.

.. admonition:: How does Open-TYNDP relate to PyPSA-Eur?

   Open-TYNDP is a soft-fork of OET/PyPSA-Eur and contains the entire Open-TYNDP project supported by OET, including code and documentation. The workflow automatically downloads publicly available data from an archived repository. OET/PyPSA-Eur itself is a soft-fork of PyPSA/PyPSA-Eur and builds on the open-source ecosystem of PyPSA.

.. admonition:: Is Open-TYNDP ready for production use?

   Open-TYNDP is under active development and is not yet feature-complete. The current development status and general limitations are important to understand before using the model. Please refer to :doc:`limitations` and :ref:`development-status` for more details. The latest version of the model is always available on `GitHub <https://github.com/open-energy-transition/open-tyndp>`__.

.. admonition:: When will Open-TYNDP be ready?

   The project is currently back-casting the 2024 TYNDP cycle to build confidence before aligning with the 2026 TYNDP cycle in Q2 2026. The :ref:`development-status` page provides a detailed roadmap of implemented and planned features.

.. admonition:: How to use Open-TYNDP?

    To explore, run, and modify Open-TYNDP, we recommend cloning the repository from `GitHub <https://github.com/open-energy-transition/open-tyndp>`__, which will give you access to the latest version of the model. As Open-TYNDP relies on a set of Python packages, you also need to install these dependencies. Please refer to :ref:`clone_repo` and :ref:`deps` for more details. However, if you are only interested in getting a first hands-on experience, we also provide a lightweight web-based experience using the `interactive workshops notebooks <https://open-energy-transition.github.io/open-tyndp-workshops>`_.


Technical Questions
===================

.. admonition:: Which operating systems are supported?

   The Open-TYNDP workflow is continuously tested for Linux and MacOS.

.. admonition:: I'm having trouble installing Open-TYNDP or getting started. Where should I start?

   The most common installation issues involve Python environment setup and solver configuration. We recommend using ``pixi`` for environment management. For solver setup, HiGHS is included by default for testing, but commercial solvers are supported as well. See :doc:`installation` for detailed platform-specific instructions, solver configuration guidance and alternative environment manager if you prefer using ``conda``.


.. admonition:: What computational resources do I need to run Open-TYNDP models?

   Full TYNDP scenario runs require significant computational resources: typically 55GB RAM, 8 CPU cores, and 1h15 runtime for NT scenario and a single planning horizon, using a commercial solver such as Gurobi. However, CBA assessment requirements are lower, typically running on standard workstations and HiGHS for around a minute per project. For testing and exploration, you can use smaller configurations using reduced temporal/spatial resolution that run on standard workstations and HiGHS. The TYNDP test configuration defined by ``config/test/config.tyndp.yaml`` is a good starting point. You can also explore lightweight and web-based examples using the `interactive workshops notebook <https://open-energy-transition.github.io/open-tyndp-workshops>`_.

.. admonition:: What solver do I need to solve Open-TYNDP models?

    It depends on the model you want to run. We recommend using HiGHS for exploring and testing the models at low temporal resolution, typically `52SEG`. HiGHS can also be used for CBA assessments. However, with higher temporal resolution, the SB models are larger and require a commercial solver.


.. admonition:: My workflow is failing or producing unexpected results. How do I troubleshoot?

   Start by running ``snakemake -call --configfile config/config.tyndp.yaml -n`` (dry-run) to validate workflow structure without execution. Then, check log files in ``logs/`` and verify intermediate results at each workflow stage. For persistent issues, see :doc:`support` for community assistance channels.

.. admonition:: I would like to develop my own features in a fork. What are your recommendations?

   We recommend following the `soft-fork strategy <https://open-energy-transition.github.io/handbook/docs/Engineering/SoftForkStrategy>`_ maintained by Open Energy Transition (OET). This approach allows you to maintain your own fork with custom features while staying synchronized with upstream improvements from both Open-TYNDP and PyPSA-Eur. The strategy provides guidance on organizing changes, managing merge conflicts, and contributing improvements back to the upstream repositories when appropriate.

Model Data and Assumptions
===========================

.. admonition:: What data sources and assumptions does Open-TYNDP use, and where can I find them?

   Open-TYNDP integrates TYNDP 2024 data (electricity demand, hydrogen topology, PECD renewable profiles, PEMMDB capacities, CBA projects and more) with PyPSA-Eur's open-source workflow. All data sources and their licenses are documented in :doc:`tyndp-2024-bundle`. The Open-TYNDP reproduces the 2024 TYNDP-specific assumptions and methodology based on the `TYNDP 2024 Scenarios Methodology Report <https://2024.entsos-tyndp-scenarios.eu/wp-content/uploads/2025/01/TYNDP_2024_Scenarios_Methodology_Report_Final_Version_250128.pdf>`__.

.. admonition:: How can I verify that Open-TYNDP results are reliable and accurate?

   Open-TYNDP includes a comprehensive benchmarking framework that validates model outputs against published TYNDP 2024 data using a multi-criteria approach. The framework is documented in :doc:`benchmarking` and the current results for NT scenario in the dedicated section of the `3rd workshop notebook <https://open-energy-transition.github.io/open-tyndp-workshops/20251203-workshop-pypsa-03.html#benchmark-results>`__. Be aware that Open-TYNDP is under active development (see :doc:`limitations` and :ref:`development-status`), and validation is ongoing as features are implemented.

Flexibility and Independence
=============================

.. admonition:: Can I use Open-TYNDP with my own private data?

   Yes, you can use Open-TYNDP with your own private data. The open-source nature of the codebase means you have full flexibility to integrate your proprietary or confidential data without any obligation to make it public. We covered this topic during the `3rd workshop <https://open-energy-transition.github.io/open-tyndp-workshops/20251203-workshop-pypsa-03.html#modify-assumptions>`__.

.. admonition:: Do I need to share my code modifications or developments?

   No, there is no obligation to share your code modifications or developments publicly. You have full governance over your own fork. You can keep your code and development private while still benefiting from updates and improvements in the Open-TYNDP repository. You can update your private codebase based on changes in Open-TYNDP at your own pace.

.. admonition:: Am I required to contribute my changes back to Open-TYNDP?

   No, contributions are welcome but entirely voluntary. You are free to use Open-TYNDP without any obligation to contribute back. However, contributing improvements, bug fixes, or new features helps strengthen the ecosystem and benefits the broader community.

.. admonition:: Can I collaborate with other organisations privately?

   Yes, you can collaborate with other organisations privately. There is no obligation to share data or code publicly when working with partners. Independent organisations can develop their own private repositories using the publicly available Open-TYNDP codebase, enabling private collaboration while maintaining interoperability through the shared foundation.

.. admonition:: What are the benefits of the open-source approach?

   The open-source approach provides several benefits: full independence in governance and decision-making, flexibility to keep parts of your work private, ability to request support and feature development from other actors in the ecosystem, interoperability with other organisations using the same foundation, and opportunities for co-developing features through partnerships when desired.

Explore results
================

.. admonition:: Where can I find the latest results?

   The preliminary results for the NT scenario and the corresponding benchmarking outputs are presented in the `3rd workshop notebook <https://open-energy-transition.github.io/open-tyndp-workshops/20251203-workshop-pypsa-03.html#benchmark-results>`_. However, new fixes and features are constantly integrated in the model. The latest networks and corresponding benchmarks are also available in a `ZIP archive <https://storage.googleapis.com/open-tyndp-data-store/runs/NT-1H-20251209.zip>`__. Be aware that Open-TYNDP is under active development (see :doc:`limitations` and :ref:`development-status`), and validation is ongoing as features are implemented.

.. admonition:: Is there a way to interactively visualize results?

   Yes, the `PyPSA-Explorer <https://github.com/open-energy-transition/PyPSA-Explorer>`_ provides interactive visualization capabilities for Open-TYNDP results. This tool was introduced in the `3rd workshop <https://open-energy-transition.github.io/open-tyndp-workshops/20251203-workshop-pypsa-03.html#interactive-exploration-with-pypsa-explorer>`_.

Contributing and Support
========================

.. admonition:: How can I contribute to Open-TYNDP?

   We strongly welcome contributions! You can file issues or make pull requests on `Github <https://github.com/open-energy-transition/open-tyndp>`_ or directly on the `PyPSA-Eur Upstream <https://github.com/PyPSA/PyPSA-Eur>`_. Please also refer to the :doc:`contributing` section.

.. admonition:: Where can I get help if I encounter issues?

   Please refer to the :doc:`support` page for various ways to reach out to us and the community, including Discord, mailing lists, and issue trackers.

.. admonition:: Where can I report bugs or request features?

   For bugs and feature requests, please use the `Open-TYNDP issues <https://github.com/open-energy-transition/open-tyndp/issues>`_.
