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
import os
import subprocess
import warnings

from pathlib import Path
from importlib.metadata import version

import pandas as pd

from barcodeqc import qc
from barcodeqc.config import output_dir_from_sample_name
from barcodeqc.files import BarcodeFileError, WildcardFileError
from barcodeqc.logging import setup_logging
from barcodeqc.utils import ExternalDependencyError

warnings.filterwarnings('ignore', category=FutureWarning)

logging.getLogger('matplotlib').setLevel(logging.WARNING)

EX_OK = getattr(os, "EX_OK", 0)
EX_USAGE = getattr(os, "EX_USAGE", 64)
EX_DATAERR = getattr(os, "EX_DATAERR", 65)
EX_NOINPUT = getattr(os, "EX_NOINPUT", 66)
EX_UNAVAILABLE = getattr(os, "EX_UNAVAILABLE", 69)
EX_SOFTWARE = getattr(os, "EX_SOFTWARE", 70)
EX_OSERR = getattr(os, "EX_OSERR", 71)
EX_IOERR = getattr(os, "EX_IOERR", 74)
EX_NOPERM = getattr(os, "EX_NOPERM", 77)


def _exit_code_for_exception(exc: Exception) -> int:
    if isinstance(exc, ExternalDependencyError):
        return EX_UNAVAILABLE
    if isinstance(exc, FileNotFoundError):
        return EX_NOINPUT
    if isinstance(exc, (BarcodeFileError, WildcardFileError, ValueError)):
        return EX_DATAERR
    if isinstance(exc, subprocess.CalledProcessError):
        return EX_IOERR
    if isinstance(exc, PermissionError):
        return EX_NOPERM
    if isinstance(exc, OSError):
        return EX_OSERR
    return EX_SOFTWARE


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="barcodeqc",
        description="Inital FASTQ-based QC of epigenomic DBiT-seq experiments from AtlasXomics.",
        epilog="Questions? Comments? Contact support@atlasxomics.com.",
    )
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"%(prog)s {version('barcodeqc')}"
    )

    subparsers = parser.add_subparsers(
        title="commands", dest="command", required=True
    )

    qc_parser = subparsers.add_parser(
        "qc",
        prog="barcodeqc qc",
        description="Inital FASTQ-based QC of epigenomic DBiT-seq experiments from AtlasXomics.",
        epilog="Questions? Comments? Contact support@atlasxomics.com.",
        help="Run barcode QC pipeline",
    )
    qc_parser.add_argument(
        "sample_name",
        help="Provide sample name for experiment."
    )
    qc_parser.add_argument(
        "r2_path",
        type=Path,
        help="Path to Read2 fastq file; accepts .fastq or .fastq.gz."
    )
    qc_parser.add_argument(
        "barcode_set",
        type=str,
        metavar="barcode_set",
        choices=[
            "bc50", "bc96", "fg96", "bc220", "bc220_05-OCT", "bc220_20-MAY"
        ],
        help="Barcode Set: bc50|bc96|fg96|bc220|bc220_05-OCT|bc220_20|"
    )
    qc_parser.add_argument(
        "-r",
        "--sample_reads",
        metavar="<reads>",
        type=int,
        required=False,
        default=10_000_000,
        help="Value to subsample reads to (default: %(default)s)."
    )
    qc_parser.add_argument(
        "-s",
        "--random_seed",
        metavar="<seed>",
        type=int,
        required=False,
        default=42,
        help="Seed for randomization during subsampling \
            (default: %(default)s)."
    )
    qc_parser.add_argument(
        "-t",
        "--tissue_position_file",
        metavar="<file>",
        type=Path,
        required=False,
        default=None,
        help="Optional tissue_positions_list.csv from AtlasXBrowser, mapping \
            barcodes to on/off tissue call and coordinates."
    )
    qc_parser.add_argument(
        "--dry_run",
        required=False,
        action="store_true",
        help="Print commands instead of executing them (default off)"
    )
    qc_parser.add_argument(
        "--count_raw_reads",
        required=False,
        action="store_true",
        help=(
            "Count total reads in input FASTQ for report metadata. "
            "Disabled by default because it scans the full file and can be "
            "slow on large datasets (default off)."
        ),
    )

    report_parser = subparsers.add_parser(
        "report",
        prog="barcodeqc report",
        help="Generate report from existing run files",
        description="Rerun report generation on results of a previous \
            execution of `barcodeqc qc`.",
        epilog="Questions? Comments? Contact support@atlasxomics.com.",
    )
    report_parser.add_argument(
        "sample_name",
        help="Sample name for report title/output."
    )
    report_parser.add_argument(
        "sample_dir",
        type=Path,
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
            raise FileNotFoundError(f"Sample directory not found: {sample_dir}")

        out_dir = output_dir_from_sample_name(args.sample_name)
        out_dir.mkdir(parents=True, exist_ok=True)

        figures_dir = sample_dir / "figures"
        tables_dir = sample_dir / "tables"
        if figures_dir.exists():
            figures = list(figures_dir.glob("*.html"))
            figures += list(figures_dir.glob("*.png"))
        else:
            figures = list(sample_dir.glob("*.html"))
            figures += list(sample_dir.glob("*.png"))
        summary_path = tables_dir / "qc_table.csv"
        onoff_path = tables_dir / "onoff_tissue_table.csv"
        if not summary_path.exists():
            summary_path = sample_dir / "qc_table.csv"
        if not onoff_path.exists():
            onoff_path = sample_dir / "onoff_tissue_table.csv"

        summary = pd.read_csv(summary_path) if summary_path.exists() else None
        onoff = pd.read_csv(onoff_path) if onoff_path.exists() else None

        linker_metrics = report.load_linker_metrics_from_dir(sample_dir)
        input_params = report.load_input_params_from_dir(sample_dir)

        if summary is not None:
            from barcodeqc.report import print_summary_table

            print_summary_table(summary)

        report.generate_report(
            figure_paths=figures,
            output_dir=out_dir,
            sample_name=args.sample_name,
            summary_table=summary,
            onoff_table=onoff,
            input_params=input_params,
            linker_metrics=linker_metrics,
            file_tag="bcQC",
            table_dir=tables_dir if tables_dir.exists() else None,
        )
        logger.info("Report generated in %s", out_dir)
        return EX_OK

    sample_dir = output_dir_from_sample_name(args.sample_name)
    sample_dir.mkdir(parents=True, exist_ok=True)

    if not args.dry_run:
        spatial_table = qc(
            sample_name=args.sample_name,
            r2_path=args.r2_path,
            barcode_set=args.barcode_set,
            sample_reads=args.sample_reads,
            random_seed=args.random_seed,
            tissue_position_file=args.tissue_position_file,
            count_raw_reads=args.count_raw_reads,
        )
        if not spatial_table.exists():
            raise RuntimeError(
                f"Spatial table not found after qc: {spatial_table}"
            )

    return EX_OK


def run(argv: list[str] | None = None) -> int:
    try:
        args = parse_args(argv)
    except SystemExit as exc:
        code = exc.code if isinstance(exc.code, int) else EX_SOFTWARE
        if code == 0:
            return EX_OK
        return EX_USAGE

    log_dir = None
    if args.command == "qc":
        log_dir = output_dir_from_sample_name(args.sample_name) / "logs"

    setup_logging(
        log_file="barcodeqc.log",
        stdout_level=logging.INFO,
        file_level=logging.DEBUG,
        log_dir=log_dir,
    )

    logger = logging.getLogger(__name__)
    try:
        return main(args)
    except KeyboardInterrupt:
        logger.error("Interrupted by user.")
        return 130
    except Exception as exc:
        exit_code = _exit_code_for_exception(exc)
        if exit_code == EX_SOFTWARE:
            logger.exception("Unhandled exception in CLI run")
        else:
            logger.error(str(exc))
        return exit_code


if __name__ == "__main__":
    raise SystemExit(run())
