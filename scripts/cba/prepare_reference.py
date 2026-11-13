# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Prepare the CBA reference network of a planning horizon.

Applies the capacity corrections from the SB to the CBA reference grid
(see `fix_reference_sb_to_cba.reconcile_reference_grid`) to the simplified SB network.
"""

import logging

import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def apply_reference_corrections(
    n: pypsa.Network,
    corrections_fn: str,
    investment_year: int,
) -> None:
    """
    Add the capacity corrections from the reference grid total to the CBA reference
    grid onto the link `p_nom` of the network.

    Parameters
    ----------
    n : pypsa.Network
        Network to correct, modified in place.
    corrections_fn : str
        Path to the per-horizon corrections CSV produced by `fix_reference_sb_to_cba`.
    investment_year : int
        Planning horizon the corrections apply to, used for logging only.

    Raises
    ------
    ValueError
        If a corrected link is not in the network or its corrected `p_nom` is negative.
    """
    corrections = pd.read_csv(corrections_fn, index_col="link_id")["p_nom"]

    missing = corrections.index.difference(n.links.index)
    if not missing.empty:
        raise ValueError(
            f"Reference grid corrections for links not in the network: {', '.join(missing)}"
        )

    p_nom = n.links.loc[corrections.index, "p_nom"] + corrections
    negative = p_nom.index[p_nom < 0]
    if not negative.empty:
        raise ValueError(
            f"Reference grid corrections lead to negative p_nom for links: {', '.join(negative)}"
        )

    logger.info(
        f"Reference grid corrections for planning horizon {investment_year} (MW):\n{corrections.to_string()}"
    )
    n.links.loc[corrections.index, "p_nom"] = p_nom


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_reference",
            horizon="2030",
            run="NT",
            configfiles=["config/config.tyndp.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)

    apply_reference_corrections(
        n, snakemake.input.corrections, snakemake.wildcards.horizon
    )

    n.export_to_netcdf(snakemake.output.network)
