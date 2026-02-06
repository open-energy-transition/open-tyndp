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
import socket
import sys
import webbrowser
from pathlib import Path
from threading import Timer

import pypsa
from pypsa_explorer import create_app

logger = logging.getLogger(__name__)

DEFAULT_PORT = 8050


def find_free_port(start_port=8050, max_attempts=50):
    """
    Find the first available port starting from start_port.
    """
    for port in range(start_port, start_port + max_attempts):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(("127.0.0.1", port))
                return port
        except OSError:
            continue
    raise RuntimeError(
        f"Could not find free port in range {start_port}-{start_port + max_attempts}"
    )


def import_network(fn: str):
    """
    Import PyPSA network and set default color for 'none' carrier.
    """
    n = pypsa.Network(fn)
    n.carriers.loc["none", "color"] = "#000000"
    return n


def open_browser(port):
    """
    Opens browser with chosen port.
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
        print("Running from Snakemake directly.")
    else:
        # Running from command line
        logging.basicConfig(level=logging.INFO)
        output_log = sys.argv[1]
        files = sys.argv[2:]
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

    # Find free port
    PORT = find_free_port(start_port=DEFAULT_PORT, max_attempts=50)
    logger.info(f"Using available port {PORT}")

    # Create the PyPSA-Explorer dash app
    app = create_app(networks)

    # Open browser with short delay
    logger.info("Launching PyPSA-Explorer...")
    Timer(1.5, open_browser, args=[PORT]).start()

    # Launch the App
    logger.info(f"Explorer launched in background at http://127.0.0.1:{PORT}")
    app.run(port=PORT, debug=False)
