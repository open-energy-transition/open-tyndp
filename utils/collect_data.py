# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Fill the local data cache with every dataset the workflow retrieves.

Run directly or through the pixi tasks, which already specify config/config.tyndp.yaml:

    # Scenario Building
    python utils/collect_data.py --configfile path/to/your/config.yaml [ARGS]
    pixi run collect-data [ARGS]
    # the CBA, plus the SB datasets it builds on
    python utils/collect_data.py --cba --configfile path/to/your/config.yaml [ARGS]
    pixi run collect-data-cba [ARGS]

Arguments other than the options below are passed to every Snakemake call, so `-n` lists
what would be downloaded without fetching it, and `-c`, `--scheduler`, `--executor` and
`--profile` reach Snakemake unchanged. The last occurrence of an argument wins, so `-c8`
overrides the `-call` used for the calls that download, and passing `--configfile` to the
pixi task replaces the one the pixi task specifies by default.
"""

import argparse
import re
import shlex
import subprocess
import sys

FILL_CACHE_CONFIG = "data={local_cache: {enable: true, fill: true}}"
READ_CACHE_CONFIG = "data={local_cache: {enable: true, fill: false}}"
RETRIEVE_RULE = re.compile(r"retrieve_[a-z0-9_]+")
DRY_RUN_FLAGS = {"-n", "--dry-run", "--dryrun"}


def run_snakemake(
    step: str, *args: str, capture: bool = False, verbose: bool = False
) -> str:
    """
    Run one Snakemake call, aborting the script if it fails.

    Captured output is held back while the call succeeds and written out when it fails, so
    Snakemake's errors always reach the terminal. Snakemake's stderr is never captured.

    Parameters
    ----------
    step : str
        What the call is doing, named in the progress line and in the error message.
    *args : str
        Command line arguments for Snakemake.
    capture : bool, default False
        Return what Snakemake wrote to stdout rather than letting it through.
    verbose : bool, default False
        Show the command and, when captured, Snakemake's output on success too.

    Returns
    -------
    str
        Snakemake's standard output, empty unless `capture` is set.

    Raises
    ------
    SystemExit
        If Snakemake exits with a non-zero return code.
    """
    command = [sys.executable, "-m", "snakemake", *args]
    print(f"\n{step}:", flush=True)
    if verbose:
        print(f"  {' '.join(command)}", flush=True)
    completed = subprocess.run(
        command,
        text=True,
        encoding="utf-8",
        errors="replace",
        stdout=subprocess.PIPE if capture else None,
    )
    if capture and (verbose or completed.returncode):
        sys.stdout.write(completed.stdout)
    if completed.returncode:
        raise SystemExit(f"{step} failed with exit code {completed.returncode}.")
    return completed.stdout or ""


def retrieve_rules(listing: str) -> list[str]:
    """
    Return the retrieve rules named in a Snakemake job listing.

    Parameters
    ----------
    listing : str
        Standard output of a Snakemake dry run.

    Returns
    -------
    list of str
        Rule names, deduplicated and sorted.
    """
    return sorted(
        {
            match.group()
            for line in listing.splitlines()
            if (match := RETRIEVE_RULE.match(line))
        }
    )


def parse_arguments() -> tuple[argparse.Namespace, list[str]]:
    """
    Read the command line, keeping unknown arguments for Snakemake.

    Returns
    -------
    argparse.Namespace
        The parsed `cba` and `configfile` options.
    list of str
        The remaining arguments, to be passed on to Snakemake.
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--cba",
        action="store_true",
        help="also collect the datasets the CBA workflow needs on top of Scenario Building",
    )
    parser.add_argument(
        "--configfile",
        nargs="+",
        default=[],
        help="config file(s) to collect for, overriding the one the pixi task passes",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="show every Snakemake call in full, including the dry runs that only list "
        "datasets; taken by this script rather than passed on to Snakemake",
    )
    return parser.parse_known_args()


def check_cba_coverage(*args: str, verbose: bool = False) -> None:
    """
    Confirm the cache holds every dataset the expanded CBA graph asks for.

    The dry run reads from the cache, as a CBA run against it does, so datasets that are
    missing and cached outputs that Snakemake would still re-retrieve are both reported.

    Parameters
    ----------
    *args : str
        Command line arguments for Snakemake.
    verbose : bool, default False
        Show the dry run behind the check.

    Raises
    ------
    SystemExit
        If the CBA workflow would still retrieve a dataset.
    """
    listing = run_snakemake(
        "Checking that the cache covers the CBA workflow in read mode",
        "cba",
        "-n",
        *args,
        capture=True,
        verbose=verbose,
    )
    if missing := retrieve_rules(listing):
        quote = subprocess.list2cmdline if sys.platform == "win32" else shlex.join
        interpreter = quote([sys.executable])
        listing_command = f"{interpreter} -m snakemake -n cba {quote(list(args))}"
        fetch_command = f"{interpreter} -m snakemake -call <paths> {quote(list(args))}"
        raise SystemExit(
            "The cache does not cover the CBA workflow, these rules would still "
            "retrieve data:\n  " + "\n  ".join(missing) + "\n\n"
            "The clean_projects inputs and the Scenario Building collection are meant to cover "
            "the CBA data files between them, so the collection and the workflow could have "
            "drifted apart. To fill the cache meanwhile, read the output paths those "
            f"jobs report:\n\n  {listing_command}\n\n"
            f"and fetch them directly, by replacing <paths>:\n\n  {fetch_command}\n\n"
            "Please report the issue under: "
            "https://github.com/open-energy-transition/open-tyndp/issues\n"
        )
    else:
        print("\nNo further retrieves required, the cache covers the CBA workflow.")


def main() -> None:
    """
    Collect every dataset the workflow retrieves into the local cache.

    Raises
    ------
    SystemExit
        If a Snakemake call fails, if the forced dry run names no retrieve rule, or if the
        cache does not cover the CBA workflow.
    """
    args, forwarded = parse_arguments()
    configfiles = ["--configfile", *args.configfile] if args.configfile else []
    base = [*configfiles, "--config", FILL_CACHE_CONFIG, *forwarded]
    read_base = [*configfiles, "--config", READ_CACHE_CONFIG, *forwarded]
    dry_run = bool(DRY_RUN_FLAGS.intersection(forwarded))

    print(
        "Filling the local cache. Working out which datasets this workflow needs requires "
        "building its full dependency graph, which can take some time before the first "
        "download starts."
    )

    if args.cba:
        # collect the CBA data (optionally incl. pre-solved networks) and run the clean_projects checkpoint
        run_snakemake(
            "Collecting the CBA project data, the pre-solved SB networks if configured, and "
            "running the clean_projects checkpoint",
            "collect_cba_data",
            "-c",
            "all",
            *base,
            "--forcerun",
            "clean_projects",
            verbose=args.verbose,
        )

    listing = run_snakemake(
        "Listing the datasets this workflow needs",
        "-n",
        "--forceall",
        *base,
        capture=True,
        verbose=args.verbose,
    )
    if not (needed := retrieve_rules(listing)):
        raise SystemExit(
            "Found no retrieve rules in the forced dry run. Either this config retrieves no "
            "datasets at all, or its job listing no longer parses. Pass '--verbose' to inspect "
            "job listing in more detail."
        )
    run_snakemake(
        f"{len(needed)} dataset(s) needed. Collecting missing ones",
        "-c",
        "all",
        *base,
        "--until",
        *needed,
        verbose=args.verbose,
    )

    if not args.cba:
        return
    if dry_run:
        print(
            "\nSkipping the coverage check, the clean_projects checkpoint has not run under "
            "a dry run and the CBA graph is not expanded yet."
        )
        return
    check_cba_coverage(*read_base, verbose=args.verbose)


if __name__ == "__main__":
    main()
