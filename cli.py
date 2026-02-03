#!/usr/bin/env python3.10
"""Local CLI wrapper for the Barcode QC workflow."""
from __future__ import annotations

from barcodeqc.cli import run


if __name__ == "__main__":
    raise SystemExit(run())
