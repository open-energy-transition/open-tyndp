# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

if config["electricity"]["base_network"] == "osm-raw":

    rule prepare_osm_network_release:
        params:
            line_types=config["lines"]["types"],
        input:
            base_network=resources("networks/base.nc"),
            stations_polygon=resources("osm-raw/build/geojson/stations_polygon.geojson"),
            buses_polygon=resources("osm-raw/build/geojson/buses_polygon.geojson"),
        output:
            buses=resources("osm-raw/release/buses.csv"),
            converters=resources("osm-raw/release/converters.csv"),
            lines=resources("osm-raw/release/lines.csv"),
            links=resources("osm-raw/release/links.csv"),
            transformers=resources("osm-raw/release/transformers.csv"),
            map=resources("osm-raw/release/map.html"),
        log:
            logs("prepare_osm_network_release.log"),
        benchmark:
            benchmarks("prepare_osm_network_release")
        threads: 1
        resources:
            mem_mb=1000,
        conda:
            "../envs/environment.yaml"
        script:
            "../scripts/prepare_osm_network_release.py"


if not config["electricity"]["pecd_renewable_profiles"]["pre_built"]["retrieve"]:

    def pecd_version(w):
        version = config_provider("electricity", "pecd_renewable_profiles", "version")(
            w
        )
        return {"pecd_raw": f"data/tyndp_2024_bundle/PECD/PECD_{version}"}

    rule prepare_pecd_release:
        params:
            cyears=config_provider(
                "electricity", "pecd_renewable_profiles", "pre_built", "cyears"
            ),
            available_pyears=config_provider(
                "electricity", "pecd_renewable_profiles", "available_years"
            ),
        input:
            unpack(pecd_version),
        output:
            pecd_prebuilt=directory(
                "data/tyndp_2024_bundle/PECD/PECD_{pecd_prebuilt_version}"
            ),
        log:
            logs("prepare_pecd_release_{pecd_prebuilt_version}.log"),
        benchmark:
            benchmarks("prepare_pecd_release_{pecd_prebuilt_version}")
        threads: 4
        resources:
            mem_mb=1000,
        script:
            "../scripts/prepare_pecd_release.py"
