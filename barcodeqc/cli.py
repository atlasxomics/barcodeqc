#!/usr/bin/env python3
"""Local CLI wrapper for the Barcode QC workflow.

This file provides a minimal local command-line entrypoint that mirrors the
parameters used by the Latch workflow. It invokes the scripts in
`./barcodeqc/` directly using subprocess. It supports a `--dry-run` flag to
print commands instead of executing them.
"""
from __future__ import annotations

import argparse
import logging
import warnings

from pathlib import Path

import pandas as pd

from barcodeqc import qc
from barcodeqc.logging import setup_logging

warnings.filterwarnings('ignore', category=FutureWarning)

logging.getLogger('matplotlib').setLevel(logging.WARNING)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="barcodeqc", prog="barcodeqc")
    subparsers = parser.add_subparsers(dest="command", required=True)

    qc_parser = subparsers.add_parser(
        "qc", help="Run barcode QC pipeline", prog="barcodeqc qc"
    )
    qc_parser.add_argument(
        "-n",
        "--sample_name",
        required=True,
        help="Provide sample name for experiment."
    )
    qc_parser.add_argument(
        "-f",
        "--r2_path",
        type=Path,
        required=True,
        help="Path to Read2 fastq file; accepts .fastq or .fastq.gz."
    )
    qc_parser.add_argument(
        "-b",
        "--barcode_set",
        type=str,
        choices=[
            "bc50", "bc96", "fg96", "bc220", "bc220_05-OCT", "bc220_20-MAY"
        ],
        required=True,
        help="Barcode Set: bc50|bc96|fg96|bc220|bc220_05-OCT|bc220_20|"
    )
    qc_parser.add_argument(
        "-r",
        "--sample_reads",
        type=int,
        required=False,
        default=10_000_000,
        help="Value to subsample reads to; default=10e6."
    )
    qc_parser.add_argument(
        "-s",
        "--random_seed",
        type=int,
        required=False,
        default=42,
        help="Seed for randomization during subsampling."
    )
    qc_parser.add_argument(
        "-t",
        "--tissue_position_file",
        type=Path,
        required=False,
        default=None,
        help="Standard tissue_positions_list.csv from AtlasXBrowser, mapping \
            barcodes to on/off tissue call and coordinates."
    )
    qc_parser.add_argument(
        "--dry_run",
        required=False,
        action="store_true",
        help="Print commands instead of executing them"
    )

    report_parser = subparsers.add_parser(
        "report",
        help="Generate report from existing run files.",
        prog="barcodeqc report",
    )
    report_parser.add_argument(
        "-n",
        "--sample_name",
        required=True,
        help="Sample name for report title/output."
    )
    report_parser.add_argument(
        "-d",
        "--sample_dir",
        type=Path,
        required=True,
        help="Directory containing existing run files (png/csv).",
    )

    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = build_parser()
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> int:
    logger = logging.getLogger(__name__)
    logger.info("CLI started")
    logger.debug("args: %s", vars(args))

    if args.command == "report":
        import barcodeqc.report as report

        sample_dir = args.sample_dir
        if not sample_dir.exists():
            logger.error("Sample directory not found: %s", sample_dir)
            return 1

        out_dir = Path.cwd() / args.sample_name
        out_dir.mkdir(parents=True, exist_ok=True)

        figures_dir = sample_dir / "figures"
        tables_dir = sample_dir / "tables"
        figures = (
            list(figures_dir.glob("*.png"))
            if figures_dir.exists()
            else list(sample_dir.glob("*.png"))
        )
        summary_path = tables_dir / "qc_table.csv"
        onoff_path = tables_dir / "onoff_tissue_table.csv"
        if not summary_path.exists():
            summary_path = sample_dir / "qc_table.csv"
        if not onoff_path.exists():
            onoff_path = sample_dir / "onoff_tissue_table.csv"

        summary = pd.read_csv(summary_path) if summary_path.exists() else None
        onoff = pd.read_csv(onoff_path) if onoff_path.exists() else None

        linker_metrics = report.load_linker_metrics_from_dir(sample_dir)

        if summary is not None:
            from barcodeqc.report import print_summary_table

            print_summary_table(summary)

        report.generate_report(
            figure_paths=figures,
            output_dir=out_dir,
            sample_name=args.sample_name,
            summary_table=summary,
            onoff_table=onoff,
            linker_metrics=linker_metrics,
            file_tag="bcQC",
            table_dir=tables_dir if tables_dir.exists() else None,
        )
        logger.info("Report generated in %s", out_dir)
        return 0

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
            logging.error(
                "Spatial table not found at %s after qc", spatial_table
            )
            return 1

    return 0


def run(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    log_dir = None
    if args.command == "qc":
        log_dir = Path.cwd() / args.sample_name / "logs"

    setup_logging(
        log_file="barcodeqc.log",
        stdout_level=logging.INFO,
        file_level=logging.DEBUG,
        log_dir=log_dir,
    )

    return main(args)


if __name__ == "__main__":
    raise SystemExit(run())
