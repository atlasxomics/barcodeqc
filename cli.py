#!/usr/bin/env python3
"""Local CLI wrapper for the Barcode QC workflow.

This file provides a minimal local command-line entrypoint that mirrors the
parameters used by the Latch workflow. It invokes the scripts in
`./barcodeqc/` directly using subprocess. It supports a `--dry-run` flag to
print commands instead of executing them.
"""
import argparse
import logging

from pathlib import Path

from barcodeqc import qc

logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="barcodeqc")

    parser.add_argument(
        "-n",
        "--sample_name",
        required=True,
        help="Provide sample name for experiment."
    )
    parser.add_argument(
        "-f",
        "--r2_path",
        type=Path,
        required=True,
        help="Path to Read2 fastq file; accepts .fastq or .fastq.gz."
    )
    parser.add_argument(
        "-b",
        "--barcode_set",
        type=str,
        choices=[
            "bc50", "bc96", "fg96", "bc220", "bc220_05-OCT", "bc220_20-MAY"
        ],
        required=True,
        help="Barcode Set: bc50|bc96|fg96|bc220|bc220_05-OCT|bc220_20|"
    )
    parser.add_argument(
        "-r",
        "--sample_reads",
        type=int,
        required=False,
        default=10_000_000,
        help="Value to subsample reads to; default=10e6."
    )
    parser.add_argument(
        "-s",
        "--random_seed",
        type=int,
        required=False,
        default=42,
        help="Seed for randomization during subsampling."
    )
    parser.add_argument(
        "-t",
        "--tissue_position_file",
        type=Path,
        required=False,
        default=None,
        help="Standard tissue_positions_list.csv from AtlasXBrowser, mapping \
            barcodes to on/off tissue call and coordinates."
    )
    parser.add_argument(
        "--dry_run",
        required=False,
        action="store_true",
        help="Print commands instead of executing them"
    )

    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = build_parser()
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> int:
    logging.info("args: %s", vars(args))

    sample_dir = Path.cwd() / args.sample_name
    sample_dir.mkdir(parents=True, exist_ok=True)

    if not args.dry_run:
        spatial_table = qc(
            sample_name=args.sample_name,
            r2_path=args.r2_path,
            barcode_set=args.barcode_set,
            sample_reads=args.sample_reads,
            random_seed=args.random_seed,
            tissue_position_file=args.tissue_position_file
        )
        if not spatial_table.exists():
            logging.warning(
                "Spatial table not found at %s after qc", spatial_table
            )

    return 0


def run(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    return main(args)


if __name__ == "__main__":
    raise SystemExit(run())
