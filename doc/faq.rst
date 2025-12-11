.. SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
..
.. SPDX-License-Identifier: CC-BY-4.0

##########################################
FAQ
##########################################

This page contains frequently asked questions about the Open-TYNDP project:

General Questions
==================

**Q: What is Open-TYNDP?**

A: Open-TYNDP is an open-source research and innovation project, which is a collaboration between `Open Energy Transition (OET) <https://openenergytransition.org/>`__ and ENTSO-E. The project aims to explore the option of a complementary open-source tool in the Ten-Year Network Development Plan (TYNDP) by building a workflow based on PyPSA-Eur. It provides one streamlined tool for the Scenario Building and the Cost-Benefit Analysis (CBA) of the TYNDP.

**Q: How does Open-TYNDP relate to PyPSA-Eur?**

A: Open-TYNDP is a soft-fork of OET/PyPSA-Eur and contains the entire Open-TYNDP project supported by OET, including code and documentation. The workflow automatically downloads publicly available data from an archived repository. OET/PyPSA-Eur itself is a soft-fork of PyPSA/PyPSA-Eur and builds on the open-source ecosystem of PyPSA.

**Q: Is Open-TYNDP ready for production use?**

A: Open-TYNDP is under active development and is not yet feature-complete. The current development status and general limitations are important to understand before using the model. Please refer to :doc:`limitations` and :ref:`development-status` for more details.

**Q: When will Open-TYNDP be ready?**

A: The project is currently back-casting the 2024 TYNDP cycle to build confidence before aligning with the 2026 TYNDP cycle in Q2 2026. The :ref:`development-status` page provides a detailed roadmap of implemented and planned features.

Technical Questions
===================

**Q: Which operating systems are supported?**

A: The Open-TYNDP workflow is continuously tested for Linux and MacOS.

**Q: I'm having trouble installing Open-TYNDP or getting started. Where should I start?**

A: The most common installation issues involve Python environment setup and solver configuration. We recommend using ``pixi`` for environment management. For solver setup, HiGHS is included by default for testing, but commercial solvers are supported as well. See :doc:`installation` for detailed platform-specific instructions and solver configuration guidance.

**Q: What computational resources do I need to run Open-TYNDP models?**

A: Full TYNDP scenario runs require significant computational resources: typically 55GB RAM, 8 CPU cores, and 1h15 runtime for NT scenario and a single planning horizon, using a commercial solver such as Gurobi. However, CBA assessment requirements are lower, typically running on standard workstations and HiGHS for around a minute per project. For testing and exploration, you can use smaller configurations using reduced temporal/spatial resolution that run on standard workstations and HiGHS. The TYNDP test configuration defined by ``config/test/config.tyndp.yaml`` is a good starting point. You can also explore lightweight and web-based examples using the `interactive workshops notebook <https://open-energy-transition.github.io/open-tyndp-workshops>`_.

**Q: My workflow is failing or producing unexpected results. How do I troubleshoot?**

A: Start by running ``snakemake -n`` (dry-run) to validate workflow structure without execution. Then, check log files in ``logs/`` and verify intermediate results at each workflow stage. For persistent issues, see :doc:`support` for community assistance channels.

**Q: I would like to develop my own features in a fork. What are your recommendations?**

A: We recommend following the `soft-fork strategy <https://open-energy-transition.github.io/handbook/docs/Engineering/SoftForkStrategy>`_ maintained by Open Energy Transition (OET). This approach allows you to maintain your own fork with custom features while staying synchronized with upstream improvements from both Open-TYNDP and PyPSA-Eur. The strategy provides guidance on organizing changes, managing merge conflicts, and contributing improvements back to the upstream repositories when appropriate.

Model Data and Assumptions
===========================

**Q: What data sources and assumptions does Open-TYNDP use, and where can I find them?**

A: Open-TYNDP integrates TYNDP 2024 data (electricity demand, hydrogen topology, PECD renewable profiles, PEMMDB capacities, CBA projects and more) with PyPSA-Eur's open-source workflow. All data sources and their licenses are documented in :doc:`tyndp-2024-bundle`. The Open-TYNDP reproduces the 2024 TYNDP-specific assumptions and methodology based on the `TYNDP 2024 Scenarios Methodology Report <https://2024.entsos-tyndp-scenarios.eu/wp-content/uploads/2025/01/TYNDP_2024_Scenarios_Methodology_Report_Final_Version_250128.pdf>`__.

**Q: How can I verify that Open-TYNDP results are reliable and accurate?**

A: Open-TYNDP includes a comprehensive benchmarking framework that validates model outputs against published TYNDP 2024 data using a multi-criteria approach. The framework is documented in :doc:`benchmarking` and the current results for NT scenario in the dedicated section of the `3rd workshop notebook <https://open-energy-transition.github.io/open-tyndp-workshops/20251203-workshop-pypsa-03.html#benchmark-results>`__. Be aware that Open-TYNDP is under active development (see :doc:`limitations` and :ref:`development-status`), and validation is ongoing as features are implemented.

Using Open-TYNDP: Flexibility and Independence
===============================================

**Q: Can I use Open-TYNDP with my own private data?**

A: Yes, you can use Open-TYNDP with your own private data. The open-source nature of the codebase means you have full flexibility to integrate your proprietary or confidential data without any obligation to make it public. We covered this topic during the `3rd workshop <https://open-energy-transition.github.io/open-tyndp-workshops/20251203-workshop-pypsa-03.html#modify-assumptions>`__.

**Q: Do I need to share my code modifications or developments?**

A: No, there is no obligation to share your code modifications or developments publicly. You can keep your code and development private while still benefiting from updates and improvements in the Open-TYNDP repository. You can update your private codebase based on changes in Open-TYNDP at your own pace.

**Q: Am I required to contribute my changes back to Open-TYNDP?**

A: No, contributions are welcome but entirely voluntary. You are free to use Open-TYNDP without any obligation to contribute back. However, contributing improvements, bug fixes, or new features helps strengthen the ecosystem and benefits the broader community.

**Q: Can I collaborate with other organisations privately?**

A: Yes, you can collaborate with other organisations privately. There is no obligation to share data or code publicly when working with partners. Independent organisations can develop their own private repositories using the publicly available Open-TYNDP codebase, enabling private collaboration while maintaining interoperability through the shared foundation.

**Q: What are the benefits of the open-source approach?**

A: The open-source approach provides several benefits: full independence in governance and decision-making, flexibility to keep parts of your work private, ability to request support and feature development from other actors in the ecosystem, interoperability with other organisations using the same foundation, and opportunities for co-developing features through partnerships when desired.

Explore results
================

**Q: Where can I find the latest results?**

A: The preliminary results for the NT scenario and the corresponding benchmarking outputs are presented in the `3rd workshop notebook <https://open-energy-transition.github.io/open-tyndp-workshops/20251203-workshop-pypsa-03.html#benchmark-results>`_. However, new fixes and features are constantly integrated in the model. The latest networks and corresponding benchmarks are also available in a `ZIP archive <https://storage.googleapis.com/open-tyndp-data-store/runs/NT-1H-20251209.zip>`__. Be aware that Open-TYNDP is under active development (see :doc:`limitations` and :ref:`development-status`), and validation is ongoing as features are implemented.

**Q: Is there a way to interactively visualize results?**

A: Yes, the `PyPSA-Explorer <https://github.com/open-energy-transition/PyPSA-Explorer>`_ provides interactive visualization capabilities for Open-TYNDP results. This tool was introduced in the `3rd workshop <https://open-energy-transition.github.io/open-tyndp-workshops/20251203-workshop-pypsa-03.html#interactive-exploration-with-pypsa-explorer>`_.

Contributing and Support
========================

**Q: How can I contribute to Open-TYNDP?**

A: We strongly welcome contributions! You can file issues or make pull requests on `Github <https://github.com/open-energy-transition/open-tyndp>`_ or directly on the `PyPSA-Eur Upstream <https://github.com/PyPSA/PyPSA-Eur>`_. Please also refer to the :doc:`contributing` section.

**Q: Where can I get help if I encounter issues?**

A: Please refer to the :doc:`support` page for various ways to reach out to the community, including Discord, mailing lists, and issue trackers.

**Q: Where can I report bugs or request features?**

A: For bugs and feature requests, please use the `Open-TYNDP issues <https://github.com/open-energy-transition/open-tyndp/issues>`_.
