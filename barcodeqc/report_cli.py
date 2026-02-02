from __future__ import annotations

import argparse
import logging

from pathlib import Path

import pandas as pd

from barcodeqc import report

logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="barcodeqc report")

    parser.add_argument(
        "-n",
        "--sample_name",
        required=True,
        help="Sample name for report title/output.",
    )
    parser.add_argument(
        "-d",
        "--sample_dir",
        type=Path,
        required=True,
        help="Directory containing existing run files (png/csv).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    sample_dir = args.sample_dir
    if not sample_dir.exists():
        logger.error("Sample directory not found: %s", sample_dir)
        return 1

    figures = list(sample_dir.glob("*.png"))
    summary_path = sample_dir / "qc_table.csv"
    onoff_path = sample_dir / "onoff_tissue_table.csv"

    summary = pd.read_csv(summary_path) if summary_path.exists() else None
    onoff = pd.read_csv(onoff_path) if onoff_path.exists() else None

    report.generate_report(
        figure_paths=figures,
        output_dir=sample_dir,
        sample_name=args.sample_name,
        summary_table=summary,
        onoff_table=onoff,
        linker_metrics=None,
        file_tag="bcQC",
    )

    logger.info("Report generated in %s", sample_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
