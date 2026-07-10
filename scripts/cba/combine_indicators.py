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


def combine_indicators_csv(input_files: list[str], output_file: str) -> None:
    """
    Concatenate multiple CSV files into one using the csv module.

    Reads the header from the first file, writes all rows from all files to the
    output, and ensures all files have the same header structure.

    Parameters
    ----------
    input_files : list[str]
        List of paths to input CSV files.
    output_file : str
        Path to output CSV file.
    """
    if not input_files:
        logger.warning("No input files provided")
        # Create empty output file with no header
        with open(output_file, "w", newline="") as f:
            pass
        return

    logger.info(f"Combining {len(input_files)} indicator files")

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

    logger.info(f"Combined {row_count} rows from {len(input_files)} files")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("combine_indicators")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Combine all indicators into a single CSV
    combine_indicators_csv(snakemake.input.indicators, snakemake.output.indicators)
