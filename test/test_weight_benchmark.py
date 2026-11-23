# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

import numpy as np
import pandas as pd
import pytest

from scripts.sb.weight_benchmark import weight_benchmarks

COLUMNS = ["table", "bus", "border", "unit", "value"]
WEIGHTS = [0.75, 0.25]


def _benchmark(values: dict[str, float | None]) -> pd.DataFrame:
    rows = [
        (table, "DE00", np.nan, "MWh" if "price" not in table else "EUR/MWh", value)
        for table, value in values.items()
        if value is not None
    ]
    return pd.DataFrame(rows, columns=COLUMNS).astype({"value": float})


@pytest.mark.parametrize(
    ("values", "expected"),
    [
        ((100.0, 200.0), 125.0),
        ((100.0, np.nan), np.nan),
        ((np.nan, 200.0), np.nan),
        ((100.0, None), np.nan),
        ((None, 200.0), np.nan),
    ],
)
def test_weight_benchmarks(values, expected):
    weighted = weight_benchmarks(
        [_benchmark({"power_generation": v}) for v in values], WEIGHTS
    )

    assert weighted.columns.tolist() == COLUMNS
    assert weighted.table.tolist() == ["power_generation"]
    assert weighted.value.item() == pytest.approx(expected, nan_ok=True)


@pytest.mark.parametrize(
    ("prices", "loads", "expected"),
    [
        ((100.0, 200.0), (10.0, 30.0), (125.0, 150.0)),
        ((100.0, np.nan), (10.0, 0.0), (np.nan, np.nan)),
        ((100.0, None), (10.0, 30.0), (np.nan, np.nan)),
        ((100.0, 200.0), (10.0, None), (125.0, np.nan)),
    ],
)
def test_weight_benchmarks_prices(prices, loads, expected):
    benchmarks = [
        _benchmark({"electricity_price": p, "electricity_price_load": load})
        for p, load in zip(prices, loads)
    ]

    weighted = weight_benchmarks(benchmarks, WEIGHTS).set_index("table").value

    assert weighted.drop("electricity_price_load").to_dict() == pytest.approx(
        {
            "electricity_price_cluster_weighted": expected[0],
            "electricity_price_load_weighted": expected[1],
        },
        nan_ok=True,
    )
