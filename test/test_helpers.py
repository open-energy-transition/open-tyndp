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
    set_scenario_config,
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


def _weather_cluster_config(tmp_path, enable=True, scenarios=True) -> dict:
    clusters = tmp_path / "clusters.csv"
    clusters.write_text(
        "target_year,cluster,representative_ws,weight,cluster_size\n"
        "2030,0,29,0.3,9\n2030,1,3,0.7,21\n2040,0,65,1.0,30\n2050,0,91,1.0,30\n"
    )
    scenario_file = tmp_path / "scenarios.yaml"
    scenario_file.write_text("NT:\n  tyndp_scenario: NT\n")
    return {
        "run": {
            "name": "NT",
            "scenarios": {"enable": scenarios, "file": scenario_file},
        },
        "planning_horizons": [2030, 2040],
        "weather_clusters": {"enable": enable, "file": clusters},
    }


def test_expand_weather_clusters(tmp_path):
    config = _weather_cluster_config(tmp_path)
    scenarios = {"NT": {"tyndp_scenario": "NT"}}

    run_parents = expand_weather_clusters(config, scenarios)

    assert run_parents == {
        "NT-py2030-ws029": "NT",
        "NT-py2030-ws003": "NT",
        "NT-py2040-ws065": "NT",
    }
    assert config["run"]["name"] == list(run_parents)
    assert scenarios["NT-py2040-ws065"] == {
        "tyndp_scenario": "NT",
        "planning_horizons": [2040],
        "wscenarios_tyndp": {2040: [65]},
    }


def test_expand_weather_clusters_disabled(tmp_path):
    config = _weather_cluster_config(tmp_path, enable=False)

    assert expand_weather_clusters(config, {"NT": {}}) == {}
    assert config["run"]["name"] == "NT"


@pytest.mark.parametrize(
    ("scenarios", "override", "match"),
    [
        (False, {}, "requires `run: scenarios: enable: true`"),
        (
            True,
            {"planning_horizons": [2030, 2035]},
            "Planning horizon 2035 of scenario 'NT' has no weather clusters",
        ),
    ],
)
def test_expand_weather_clusters_raises(tmp_path, scenarios, override, match):
    config = _weather_cluster_config(tmp_path, scenarios=scenarios)

    with pytest.raises(ValueError, match=match):
        expand_weather_clusters(config, {"NT": override})


@pytest.mark.parametrize(
    ("enable", "run", "horizons", "wscenarios"),
    [
        (True, "NT-py2030-ws029", [2030], {2030: [29]}),
        (True, "NT", [2030, 2040], None),
        (False, "NT", [2030, 2040], None),
    ],
)
def test_set_scenario_config(tmp_path, enable, run, horizons, wscenarios):
    config = _weather_cluster_config(tmp_path, enable=enable)

    set_scenario_config(
        SimpleNamespace(config=config, wildcards=Wildcards(fromdict={"run": run}))
    )

    assert config["tyndp_scenario"] == "NT"
    assert config["planning_horizons"] == horizons
    assert config.get("wscenarios_tyndp") == wscenarios


def test_set_scenario_config_disabled_child(tmp_path):
    config = _weather_cluster_config(tmp_path, enable=False)
    wildcards = Wildcards(fromdict={"run": "NT-py2030-ws029"})

    with pytest.raises(KeyError):
        set_scenario_config(SimpleNamespace(config=config, wildcards=wildcards))
