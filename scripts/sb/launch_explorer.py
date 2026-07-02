# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT
"""
Launches the PyPSA-Explorer web interface to visualize workflow result networks.
The explorer will automatically use the next available port starting from 8050.

To close all running explorer instances when finished, run:
`snakemake -c1 close_explorers --configfile config/config.tyndp.yaml`
"""

import logging
import re
import sys
import webbrowser
from pathlib import Path
from threading import Timer

import pypsa
from pypsa_explorer import create_app

logger = logging.getLogger(__name__)

# Allow-list of characters permitted in a command-line file path
ALLOWED_PATH = re.compile(r"[A-Za-z0-9._/\\-]+")

# Highest valid TCP port (2**16 - 1)
MAX_PORT = 65535


def sanitize_path(path_str: str, label: str) -> Path:
    """
    Validate an untrusted path and resolve it inside the working directory.

    Guards against path manipulation: only allow-listed characters are accepted,
    absolute paths and ``..`` traversal are rejected, and the resolved path must
    stay within the current working directory.

    Parameters
    ----------
    path_str : str
        Untrusted path, expected to be relative to the working directory.
    label : str
        Human-readable input name, used in error messages.

    Returns
    -------
    pathlib.Path
        The validated path resolved against the working directory.

    Raises
    ------
    ValueError
        If the path is empty, contains disallowed characters, is absolute or
        contains ``..``, or resolves outside the working directory.
    """
    if not path_str or not ALLOWED_PATH.fullmatch(path_str):
        raise ValueError(f"Invalid characters in {label}: {path_str!r}")
    path = Path(path_str)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"{label} must be relative and without '..': {path_str!r}")
    base_dir = Path.cwd().resolve()
    resolved = (base_dir / path).resolve()
    if not resolved.is_relative_to(base_dir):
        raise ValueError(f"{label} escapes the working directory: {path_str!r}")
    return resolved


def sanitize_port(port_str: str) -> int:
    """
    Parse an untrusted TCP port and ensure it lies in the valid range.

    Parameters
    ----------
    port_str : str
        Untrusted port value to validate.

    Returns
    -------
    int
        The parsed port number, guaranteed to be within 1-``MAX_PORT``.

    Raises
    ------
    ValueError
        If the value is not an integer within the range 1-``MAX_PORT``.
    """
    if not (port_str.isdigit() and 1 <= (port := int(port_str)) <= MAX_PORT):
        raise ValueError(f"Port must be an integer in 1-{MAX_PORT}, got: {port_str!r}")
    return port


def sanitize_input_files(path_strs: list[str]) -> list[str]:
    """
    Validate untrusted network file paths.

    Each path must resolve inside the working directory (see ``sanitize_path``),
    carry a ``.nc`` suffix, and refer to an existing file.

    Parameters
    ----------
    path_strs : list[str]
        Untrusted network file paths, expected to be relative to the working
        directory.

    Returns
    -------
    list[str]
        The validated paths resolved against the working directory.

    Raises
    ------
    ValueError
        If any path fails validation, lacks a ``.nc`` suffix, or does not exist.
    """
    files = []
    for path_str in path_strs:
        path = sanitize_path(path_str, "network file")
        if path.suffix != ".nc" or not path.is_file():
            raise ValueError(f"Not an existing '.nc' file: {path_str!r}")
        files.append(str(path))
    return files


def import_network(fn: str):
    """
    Import PyPSA network and set default color for 'none' carrier.

    Parameters
    ----------
    fn : str
        Path to the PyPSA network file to load.

    Returns
    -------
    pypsa.Network
        The loaded network with a default color assigned to the 'none' carrier.
    """
    n = pypsa.Network(fn)
    n.carriers.loc["none", "color"] = "#000000"
    return n


def open_browser(port):
    """
    Opens browser with chosen port.

    Parameters
    ----------
    port : int
        Port on which the PyPSA-Explorer app is launched.

    Returns
    -------
    None
    """
    webbrowser.open_new(f"http://127.0.0.1:{port}")


if __name__ == "__main__":
    # Check if running from command line (has arguments) or from Snakemake
    running_from_cmdline = len(sys.argv) > 1
    if not running_from_cmdline and "snakemake" not in globals():
        # For debugging
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "launch_explorer",
            configfiles="config/config.tyndp.yaml",
            run="NT-cy2009-20260130",
        )

    # Get files and output path from snakemake or command line arguments
    if "snakemake" in globals():
        # Running from Snakemake directly for debugging
        # Add parent directory to path to find scripts module
        sys.path.insert(0, str(Path(__file__).parent.parent.parent))
        from scripts._helpers import configure_logging

        configure_logging(snakemake)
        files = snakemake.input
        output_log = snakemake.output[0]
        port = snakemake.params.port
        print("Running from Snakemake directly.")
    else:
        # Running from command line
        logging.basicConfig(level=logging.INFO)
        output_log = sanitize_path(sys.argv[1], "log path")
        port = sanitize_port(sys.argv[2])
        files = sanitize_input_files(sys.argv[3:])
        print("Running from command line.")

    # Add file handler to write to explorer_launched.log
    file_handler = logging.FileHandler(output_log)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    )
    logger.addHandler(file_handler)

    # Load networks into a dictionary for PyPSA-Explorer
    networks = {fn.split("_")[-1].split(".")[0]: import_network(fn) for fn in files}
    logger.info(f"Successfully loaded {len(networks)} networks: {networks}")

    # Create the PyPSA-Explorer dash app
    app = create_app(networks)

    # Open browser with short delay
    logger.info("Launching PyPSA-Explorer...")
    logger.info(f"Using available port {port}.")
    Timer(1.5, open_browser, args=[port]).start()

    # Launch the App
    logger.info(f"Explorer launched in background at http://127.0.0.1:{port}")
    app.run(port=port, debug=False)
