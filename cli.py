#!/usr/bin/env python3
"""Local CLI wrapper for the Barcode QC workflow.

This file provides a minimal local command-line entrypoint that mirrors the
parameters used by the Latch workflow. It invokes the scripts in
`wf/scripts` directly using subprocess. It supports a `--dry-run` flag to
print commands instead of executing them.

Notes:
- Defaults were chosen to be safe for local runs; scripts may still require
  additional environment-specific packages. Use `--dry-run` to validate
  command construction before running.
"""
import click
import logging
import os
import shlex
import subprocess
import sys

logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")


def run_command(cmd, dry_run=False, cwd=None):
    pretty = cmd if isinstance(cmd, str) else " ".join(cmd)
    logging.info("Command: %s", pretty)
    if dry_run:
        return
    subprocess.run(
        cmd if isinstance(cmd, (list, tuple))
        else shlex.split(cmd), check=True, cwd=cwd
    )


@click.command()
@click.option("--sample_name", required=True, help="Sample name")
@click.option("--remoteReadTwo", required=True, help="Path to Read2 file")
@click.option("--bcSet", default="bc220", help="Barcode set (bc220, fg96, bc96, bc50)")
@click.option("--outReads", default=10000000, type=int)
@click.option("--seed", default=100, type=int)
@click.option("--tissuePos_file", default=None, help="Path to tissue positions file")
@click.option("--output_directory", default="./outputs", help="Root output directory")
@click.option(
    "--dry-run", is_flag=True, help="Print commands instead of executing them"
)
def main(
    sample_name,
    remoteReadTwo,
    bcSet,
    outReads,
    seed,
    tissuePos_file,
    output_directory,
    dry_run,
):
    repo_root = os.path.abspath(os.path.dirname(__file__))
    scripts_dir = os.path.join(repo_root, "barcodeqc")

    sample_dir = os.path.join(os.path.abspath(output_directory), sample_name)
    os.makedirs(sample_dir, exist_ok=True)

    merchecker = os.path.join(scripts_dir, "merchecker.py")

    cmd1 = [
        sys.executable,
        merchecker,
        "-b",
        bcSet,
        "-f",
        remoteReadTwo,
        "-r",
        str(outReads),
        "-s",
        str(seed),
        "-n",
        sample_name,
        "-t",
        tissuePos_file if tissuePos_file else "none"
    ]
    run_command(cmd1, dry_run=dry_run)

    spatial_table = os.path.join(sample_dir, f"{sample_name}_spatialTable.csv")
    if not os.path.exists(spatial_table):
        logging.warning(
            "Spatial table not found at %s after merChecker", spatial_table
        )


if __name__ == "__main__":
    main()
