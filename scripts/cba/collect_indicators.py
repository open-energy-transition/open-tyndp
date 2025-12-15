# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Concatenate indicators from individual indicators files.

This script reads all individual project indicator CSV files and combines them
into a single CSV file.
"""

import csv
import logging

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def collect_indicators_csv(input_files, output_file):
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

    # Write output file
    row_count = 0
    with open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(header)

        for input_file in input_files:
            with open(input_file, newline="") as infile:
                reader = csv.reader(infile)

                # Read and verify header
                file_header = next(reader, None)
                if file_header != header:
                    logger.warning(
                        f"Header mismatch in {input_file}. "
                        f"Expected: {header}, Got: {file_header}"
                    )

                # Write all data rows
                for row in reader:
                    writer.writerow(row)
                    row_count += 1

    logger.info(f"Collected {row_count} rows from {len(input_files)} files")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("collect_indicators")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Collect all indicators into a single CSV
    collect_indicators_csv(snakemake.input.indicators, snakemake.output.indicators)
