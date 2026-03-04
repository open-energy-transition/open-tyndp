# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Collect indicators from individual indicators files and create weighted average
indicator.

This script reads all collected project indicators from the different weather years
and calculates a new combined indiactor which is stores in a new CSV file.
"""

import csv
import logging

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

INDICATOR_UNITS = {
    "B1_total_system_cost_change": "Meuro/year",
    "B2a_co2_variation": "t/year",
    "B2a_societal_cost_variation": "Meuro/year",
    "B3a_res_capacity_change": "MW",
    "B3_res_generation_change": "MWh/year",
    "B3_annual_avoided_curtailment": "MWh/year",
    "B4a_nox": "kg/year",
    "B4b_nh3": "kg/year",
    "B4c_sox": "kg/year",
    "B4d_pm25": "kg/year",
    "B4e_pm10": "kg/year",
    "B4f_nmvoc": "kg/year",
}

# TODO read from CSV file
cyear_weightings = {
    1995: 0.233,
    2008: 0.367,
    2009: 0.400,
}


def average_indicators_csv(input_files, output_file):
    """
    Concatenate multiple CSV files into one using the csv module.

    Args:
        input_files: List of paths to input CSV files
        output_file: Path to output CSV file

    The function:
    1. Reads the header from the first file
    2. Writes all rows from all files to the output
    3. Ensures all files have the same header structure
    """
    if not input_files:
        logger.warning("No input files provided")
        # Create empty output file with no header
        with open(output_file, "w", newline="") as f:
            pass
        return

    logger.info(f"Collecting {len(input_files)} indicator files")

    # Read header from first file
    with open(input_files[0], newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        # add additional field to the new header structure
        header.insert(0, "cyear")
        header.insert(0, "cyear_weight")

    # Write output file
    row_count = 0
    with open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(header)

        for input_file in input_files:
            # check if input_file shows reference to required weather year
            cyear_match = False
            cyear = None
            cyear_weight = 1
            for cy in cyear_weightings:
                cyear_str = "cy" + str(cy)
                if cyear_str in input_file:
                    cyear_match = True
                    cyear = cy
                    cyear_weight = cyear_weightings[cy]

            if cyear_match:
                logger.info(f"Process file {input_file} ...")
                with open(input_file, newline="") as infile:
                    reader = csv.reader(infile)

                    # Read and verify header
                    file_header = next(reader, None)
                    # add additional field to the new header structure
                    file_header.insert(0, "cyear")
                    file_header.insert(0, "cyear_weight")
                    if file_header != header:
                        logger.warning(
                            f"Header mismatch in {input_file}. "
                            f"Expected: {header}, Got: {file_header}"
                        )

                    # Write all data rows
                    for row in reader:
                        if "Open-TYNDP" in row:
                            row.insert(0, str(cyear))
                            row.insert(0, cyear_weight)
                        else:
                            row.insert(0, "")  # cyear
                            row.insert(0, 1)  # cyear_weight
                        writer.writerow(row)
                        row_count += 1
            else:
                logger.info(
                    f"Problem in catching cyear and cyear_weight: cy={cy}, cyear={cyear}, cyear_weight={cyear_weight}"
                )

    logger.info(f"Collected {row_count} rows from {len(input_files)} files")

    df = pd.read_csv(output_file, index_col=False)
    row_count = 0
    project_id = df.project_id.unique()[0]
    method = df.method.unique()[0]
    alternative_subindex = {"min": "low", "mean": "central", "max": "high"}

    # loop through all coolected indicators
    for INDICATOR_UNIT in INDICATOR_UNITS:
        units = df[(df.indicator == INDICATOR_UNIT)].units.unique()[0]
        indicator_values = {
            "min": (
                df[(df.indicator == INDICATOR_UNIT) & (df.source == "Open-TYNDP")].value
            ).min(),
            "mean": (
                df[(df.indicator == INDICATOR_UNIT) & (df.source == "Open-TYNDP")].value
                * df[
                    (df.indicator == INDICATOR_UNIT) & (df.source == "Open-TYNDP")
                ].cyear_weight
            ).sum(),
            "max": (
                df[(df.indicator == INDICATOR_UNIT) & (df.source == "Open-TYNDP")].value
            ).max(),
        }

        for subindex in ["min", "mean", "max"]:
            if INDICATOR_UNIT == "B2a_societal_cost_variation":
                subindex_to_use = alternative_subindex[subindex]
            else:
                subindex_to_use = subindex

            df.loc[len(df)] = dict(
                {
                    "cyear_weight": 1.0,
                    "project_id": project_id,
                    "method": method,
                    "source": "Open-TYNDP",
                    "indicator": INDICATOR_UNIT,
                    "subindex": subindex_to_use,
                    "units": units,
                    "value": indicator_values[subindex],
                }
            )
            row_count += 1

    df.to_csv(output_file, index=False)
    logger.info(f"Added {row_count} rows with averaged indices")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("collect_indicators")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Collect all indicators into a single CSV
    average_indicators_csv(snakemake.input.indicators, snakemake.output.indicators)
