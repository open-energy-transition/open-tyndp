# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Unified script for managing SLURM node exclusions.
Can be run in two modes:
1. Return current exclusions (--exclude mode)
2. Run the continuous updater loop (--update mode)
"""

import argparse
import os
import subprocess
import time


def cmd(x):
    return subprocess.check_output(x, text=True).strip()


def is_user_in_group(username: str, group: str = "rise-oet") -> bool:
    """
    Check if a user is part of a specific group.
    """
    try:
        groups_output = subprocess.check_output(
            ["id", "-nG", username.strip()], text=True
        ).strip()
        groups = groups_output.split()
        return group in groups
    except subprocess.CalledProcessError:
        return False


def get_current_exc_nodes(jobid: str) -> set[str]:
    """
    Return the current ExcNodeList for a job as a set of hostnames.
    If no exclusion is set, returns an empty set.
    """
    # Query job description
    out = subprocess.check_output(["scontrol", "show", "job", jobid], text=True)
    # Extract ExcNodeList=... token (may be absent)
    exc = ""
    for token in out.replace("\n", " ").split():
        if token.startswith("ExcNodeList="):
            exc = token.split("=", 1)[1]
            break
    if not exc or exc in ("(null)", "None", "n/a"):
        return set()
    # Expand to full hostnames
    try:
        expanded = subprocess.check_output(
            ["scontrol", "show", "hostnames", exc], text=True
        ).split()
    except subprocess.CalledProcessError:
        # If SLURM can't expand (e.g., malformed), treat as no exclusions
        expanded = []
    return set(expanded)


def expand_nodestr_to_set(nodelist_str: str) -> set[str]:
    """
    Expand a compressed SLURM nodelist string (e.g., 'htc-cmp[001-003]') to a set of hosts.
    Returns empty set if the string is empty.
    """
    if not nodelist_str:
        return set()
    expanded = subprocess.check_output(
        ["scontrol", "show", "hostnames", nodelist_str], text=True
    ).split()
    return set(expanded)


# -----------------------------
# Get nodes running "smk-solv" jobs
# -----------------------------
def get_smk_solv_nodes():
    # Get all running jobs with user info
    out = cmd(["squeue", "-h", "-o", "%.18i %.30j %.2t %R %.20u"])

    nodelists = []
    for line in out.splitlines():
        jobid, name, state, nodelist, job_user = line.split(maxsplit=4)

        if state != "R":
            continue
        if "smk-solv" not in name:
            continue
        if not nodelist.startswith("htc"):
            continue
        # Only consider jobs from users in the oet-rise group
        if not is_user_in_group(job_user):
            continue

        nodelists.append(nodelist)

    # Expand compressed node lists
    nodes = set()
    for nl in nodelists:
        expanded = cmd(["scontrol", "show", "hostnames", nl])
        for host in expanded.splitlines():
            nodes.add(host)

    return sorted(nodes)


# -----------------------------
# Compress list of nodes back to SLURM syntax
# -----------------------------
def compress_nodelist(nodes):
    if not nodes:
        return ""

    host_arg = ",".join(nodes)
    out = cmd(["scontrol", "show", "hostlist", host_arg])
    return out.strip()


# -----------------------------
# Get queued (PD) jobs from oet-rise group members
# -----------------------------
def get_queued_jobs(user=None):
    # Get all queued jobs with user info
    if user:
        out = cmd(["squeue", "-h", "-o", "%.18i %.30j %.2t", "-u", user])
    else:
        out = cmd(["squeue", "-h", "-o", "%.18i %.30j %.2t"])

    queued = []
    for line in out.splitlines():
        jobid, name, state = line.split(maxsplit=3)

        if state != "PD":
            continue

        if not name.startswith("smk"):  # Only snakemake jobs
            continue

        queued.append((jobid, name))

    return queued


# -----------------------------
# Update queued jobs
# -----------------------------
def update_jobs(excnodes_compressed: str, user=None):
    desired_set = expand_nodestr_to_set(excnodes_compressed)

    queued = get_queued_jobs(user=user)
    if not queued:
        print("No queued jobs to consider.")
        return

    for jobid, name in queued:
        current_set = get_current_exc_nodes(jobid)

        if current_set == desired_set:
            continue

        # If desired is empty and current is empty → already handled above.
        # If desired is empty but current non-empty, clear it by setting ExcNodeList=
        exc_arg = f"ExcNodeList={excnodes_compressed}"
        print(
            f"Updating job {jobid} ({name}) → {exc_arg} (current had {len(current_set)} nodes)"
        )
        subprocess.run(["scontrol", "update", f"job={jobid}", exc_arg], check=False)


# -----------------------------
# Mode: Return exclusions
# -----------------------------
def get_exclusions():
    """
    Get nodes running smk-solv jobs and print them as comma-separated list.
    This mode is equivalent to the old exclude.py functionality.
    Only considers jobs from users in the oet-rise group.
    """
    squeue_cmd = ["squeue", "-h", "-o", "%.18i|%.30j|%.2t|%R|%.20u"]
    out = subprocess.check_output(squeue_cmd, text=True)
    nodelist = []
    for line in out.splitlines():
        jobid, name, state, node, job_user = line.split("|", 4)
        if "smk-solv" not in name or not node.startswith("htc"):
            continue
        # Only consider jobs from users in the oet-rise group
        if not is_user_in_group(job_user):
            continue
        nodelist.append(node)

    print(",".join(nodelist))


# -----------------------------
# Mode: Run updater loop
# -----------------------------
def run_updater_loop(interval=2):
    """
    Continuously update exclusions at specified interval.
    This mode is equivalent to the old update_exclusions.py main() functionality.
    """

    print(f"Starting exclusions updater (interval: {interval}s)")
    print("Press Ctrl+C to stop")
    user = os.environ.get("USER", None)
    try:
        while True:
            nodes = get_smk_solv_nodes()
            exc = compress_nodelist(nodes)

            print("\nDetected smk-solv nodes:", nodes)
            print("Compressed:", exc)

            update_jobs(exc, user=user)

            time.sleep(interval)
    except KeyboardInterrupt:
        print("\nExclusions updater stopped")


# -----------------------------
# Main entry point
# -----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Manage SLURM node exclusions - either return current exclusions or run updater loop"
    )
    parser.add_argument(
        "--exclude",
        action="store_true",
        help="Return comma-separated list of nodes to exclude (default behavior)",
    )
    parser.add_argument(
        "--update", action="store_true", help="Run the continuous updater loop"
    )
    parser.add_argument(
        "--interval",
        type=int,
        default=2,
        help="Update interval in seconds (default: 2, only used with --update)",
    )

    args = parser.parse_args()

    # Default behavior is --exclude if neither flag is specified
    if args.update:
        run_updater_loop(args.interval)
    else:
        get_exclusions()


if __name__ == "__main__":
    main()
