# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Tests for scripts/cba/clean_projects.py project filtering.
"""

import logging

import pandas as pd

from scripts.cba.clean_projects import remove_unclear_border


def test_remove_unclear_border_reports_causes_separately(caplog):
    """Unparseable borders and unknown bus codes are dropped and warned about separately."""
    projects = pd.DataFrame(
        {
            "project_id": [1, 2, 3],
            "project_name": ["cross-border", "internal", "third country"],
            "is_crossborder": [True, False, True],
            "border": ["DE00-NL00", "internalDE00", "BG00-TR00"],
            "bus0": ["DE00", None, "BG00"],
            "bus1": ["NL00", None, "TR00"],
        }
    )

    with caplog.at_level(logging.WARNING):
        kept = remove_unclear_border(projects, pd.Index(["DE00", "NL00", "BG00"]))

    assert kept["project_id"].tolist() == [1]

    unparsed, unknown_bus = (r.getMessage() for r in caplog.records)
    assert "1 out of 3 project borders that are not reported as" in unparsed
    assert "internalDE00" in unparsed
    assert (
        "1 out of 3 project borders that have bus codes that are missing" in unknown_bus
    )
    assert "BG00-TR00" in unknown_bus
