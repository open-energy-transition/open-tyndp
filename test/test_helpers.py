# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from types import SimpleNamespace

import pandas as pd
import pypsa
import pytest
from pydantic import ValidationError
from snakemake.iocontainers import Wildcards

from scripts._helpers import (
    expand_weather_clusters,
    get_temporal_resolution,
    rename_network_component,
    select_weather_clusters,
    set_scenario_config,
    weather_cluster_weights,
)
from scripts.lib.validation.config.clustering import _TemporalConfig


def _build_network_with_generators() -> pypsa.Network:
    network = pypsa.Network()
    network.add("Bus", "bus")
    network.snapshots = pd.RangeIndex(2)
    network.add("Generator", "G1", bus="bus", p_nom=10)
    network.add("Generator", "G2", bus="bus", p_nom=20)
    network.generators_t.p = pd.DataFrame(
        [[0.5, 1.0], [0.7, 1.2]],
        index=network.snapshots,
        columns=["G1", "G2"],
    )
    return network


def test_rename_network_component_updates_static_and_dynamic_tables():
    network = _build_network_with_generators()
    rename_map = pd.Series(["G1_new", "G2_new"], index=["G1", "G2"])

    rename_network_component(network, "Generator", rename_map)

    static_index = network.static("Generator").index
    assert "G1_new" in static_index and "G1" not in static_index
    assert "G2_new" in static_index and "G2" not in static_index

    generator_p = network.generators_t.p
    assert list(generator_p.columns) == ["G1_new", "G2_new"]

    dynamic_tables = network.dynamic("Generator")
    assert generator_p is dynamic_tables["p"]


def test_rename_network_component_detects_name_conflicts():
    network = _build_network_with_generators()
    rename_map = pd.Series(["G2"], index=["G1"])

    with pytest.raises(ValueError, match="collide"):
        rename_network_component(network, "Generator", rename_map)


def test_rename_network_component_rejects_missing_keys():
    network = _build_network_with_generators()
    rename_map = pd.Series(["G3_new"], index=["missing"])

    with pytest.raises(KeyError, match="missing keys"):
        rename_network_component(network, "Generator", rename_map)


@pytest.mark.parametrize(
    ("temporal", "expected"),
    [
        ({"averaging": "24h"}, ("averaging", "24h")),
        ({"averaging": "1d"}, ("averaging", "1d")),
        ({"segmentation": 4380}, ("segmentation", 4380)),
        ({"representative": 3}, ("representative", 3)),
        ({}, None),
    ],
)
def test_get_temporal_resolution(temporal, expected):
    config = _TemporalConfig(**temporal)
    assert get_temporal_resolution(config.model_dump()) == expected


@pytest.mark.parametrize(
    "temporal",
    [{"averaging": "banana"}, {"averaging": "24h", "segmentation": 4380}],
)
def test_temporal_config_rejects_invalid_settings(temporal):
    with pytest.raises(ValidationError):
        _TemporalConfig(**temporal)


@pytest.fixture
def weather_clusters(tmp_path) -> dict:
    clusters = tmp_path / "clusters.csv"
    clusters.write_text(
        "horizon,cluster,weather_scenario,weight,cluster_size\n"
        "2030,0,29,0.3,9\n2030,1,3,0.6,18\n2030,2,21,0.1,3\n2040,0,65,1.0,30\n"
    )
    return {"file": str(clusters), "select": "all"}


@pytest.fixture
def tyndp_config(tmp_path, weather_clusters) -> dict:
    scenario_file = tmp_path / "scenarios.yaml"
    scenario_file.write_text(
        "NT:\n  tyndp_scenario: NT\n"
        "HEV:\n  planning_horizons: [2040]\n"
        "  weather_clusters_tyndp:\n    select: first\n"
    )
    return {
        "run": {
            "name": "NT",
            "scenarios": {"enable": True, "file": str(scenario_file)},
        },
        "tyndp_scenario": "NT",
        "planning_horizons": [2030, 2040],
        "weather_clusters_tyndp": weather_clusters,
    }


@pytest.mark.parametrize(
    ("select", "expected"),
    [
        ("all", [(2030, 29), (2030, 3), (2030, 21), (2040, 65)]),
        ("first", [(2030, 29), (2040, 65)]),
        ({2030: [3, 29], 2040: [65]}, [(2030, 3), (2030, 29), (2040, 65)]),
    ],
)
def test_select_weather_clusters(weather_clusters, select, expected):
    selected = select_weather_clusters(
        {**weather_clusters, "select": select}, [2030, 2040]
    )

    assert list(zip(selected.horizon, selected.weather_scenario)) == expected


@pytest.mark.parametrize(
    ("select", "horizons", "match"),
    [
        ({2030: [4], 2040: [65]}, [2030, 2040], r"\[\(2030, 4\)\].* are not in"),
        ({2030: [3], 2050: [91]}, [2030], r"\[\(2050, 91\)\].* are not in"),
        ({2030: [3]}, [2030, 2040], r"for planning horizons \[2040\]"),
        ("all", [2030, 2050], r"for planning horizons \[2050\]"),
    ],
)
def test_select_weather_clusters_raises(weather_clusters, select, horizons, match):
    with pytest.raises(ValueError, match=match):
        select_weather_clusters({**weather_clusters, "select": select}, horizons)


@pytest.mark.parametrize(
    ("select", "expected"),
    [
        ("all", {29: 0.3, 3: 0.6, 21: 0.1}),
        ("first", {29: 1.0}),
        ({2030: [3, 21]}, {3: 6 / 7, 21: 1 / 7}),
    ],
)
def test_weather_cluster_weights(weather_clusters, select, expected):
    weights = weather_cluster_weights({**weather_clusters, "select": select}, 2030)

    assert weights.to_dict() == pytest.approx(expected)


@pytest.mark.parametrize(
    ("name", "runs"),
    [
        (
            "NT",
            [
                "NT-py2030-ws029",
                "NT-py2030-ws003",
                "NT-py2030-ws021",
                "NT-py2040-ws065",
            ],
        ),
        (["HEV", "NT-py2030-ws003"], ["HEV-py2040-ws065", "NT-py2030-ws003"]),
    ],
)
def test_expand_weather_clusters(tyndp_config, weather_clusters, name, runs):
    tyndp_config["run"]["name"] = name
    hev_clusters = {"file": weather_clusters["file"], "select": "first"}
    scenarios = {
        "NT": {"tyndp_scenario": "NT"},
        "HEV": {"planning_horizons": [2040], "weather_clusters_tyndp": hev_clusters},
    }

    run_parents = expand_weather_clusters(tyndp_config, scenarios)

    assert run_parents == {
        "NT-py2030-ws029": "NT",
        "NT-py2030-ws003": "NT",
        "NT-py2030-ws021": "NT",
        "NT-py2040-ws065": "NT",
        "HEV-py2040-ws065": "HEV",
    }
    assert tyndp_config["run"]["name"] == runs
    assert scenarios["NT-py2030-ws003"] == {
        "tyndp_scenario": "NT",
        "planning_horizons": [2030],
        "weather_clusters_tyndp": {"select": {2030: [3]}},
    }
    assert scenarios["HEV-py2040-ws065"]["weather_clusters_tyndp"] == {
        "file": weather_clusters["file"],
        "select": {2040: [65]},
    }


def test_expand_weather_clusters_non_tyndp(tyndp_config):
    tyndp_config["tyndp_scenario"] = False
    scenarios = {"NT": {}}

    assert expand_weather_clusters(tyndp_config, scenarios) == {}
    assert scenarios == {"NT": {}}
    assert tyndp_config["run"]["name"] == "NT"


@pytest.mark.parametrize(
    ("run", "error", "match"),
    [
        ({"scenarios": {"enable": False}}, ValueError, "requires `run: scenarios"),
        (
            {"scenarios": {"enable": True, "file": "missing.yaml"}},
            FileNotFoundError,
            "missing.yaml does not exist",
        ),
        ({"name": ["NT", "NT-py2035-ws029"]}, ValueError, r"\['NT-py2035-ws029'\]"),
    ],
)
def test_expand_weather_clusters_raises(tyndp_config, run, error, match):
    tyndp_config["run"].update(run)

    with pytest.raises(error, match=match):
        expand_weather_clusters(tyndp_config, {"NT": {"tyndp_scenario": "NT"}})


@pytest.mark.parametrize(
    ("tyndp", "run", "horizons", "select"),
    [
        ("NT", "NT-py2030-ws003", [2030], {2030: [3]}),
        ("NT", "NT", [2030, 2040], "all"),
        (False, "NT", [2030, 2040], "all"),
    ],
)
def test_set_scenario_config(tyndp_config, tyndp, run, horizons, select):
    tyndp_config["tyndp_scenario"] = tyndp

    set_scenario_config(
        SimpleNamespace(config=tyndp_config, wildcards=Wildcards(fromdict={"run": run}))
    )

    assert tyndp_config["tyndp_scenario"] == "NT"
    assert tyndp_config["planning_horizons"] == horizons
    assert tyndp_config["weather_clusters_tyndp"]["select"] == select


def test_set_scenario_config_non_tyndp_child(tyndp_config):
    tyndp_config["tyndp_scenario"] = False
    wildcards = Wildcards(fromdict={"run": "NT-py2030-ws003"})

    with pytest.raises(KeyError):
        set_scenario_config(SimpleNamespace(config=tyndp_config, wildcards=wildcards))
