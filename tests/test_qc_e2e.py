from __future__ import annotations

import os
from pathlib import Path
from shutil import which

import pytest

from barcodeqc.qc import qc


@pytest.mark.e2e
def test_qc_smoke_with_local_example_data(
    monkeypatch,
    tmp_path: Path,
    example_fastq_path: Path,
    example_tissue_positions_path: Path,
) -> None:
    if os.environ.get("BARCODEQC_RUN_E2E") != "1":
        pytest.skip("set BARCODEQC_RUN_E2E=1 to run end-to-end QC smoke test")
    if which("seqtk") is None or which("cutadapt") is None:
        pytest.skip("seqtk and cutadapt must be available on PATH")

    monkeypatch.chdir(tmp_path)

    spatial_table = qc(
        sample_name="example",
        r2_path=example_fastq_path,
        barcode_set="bc220_20-MAY",
        sample_reads=100_000,
        random_seed=42,
        tissue_position_file=example_tissue_positions_path,
        count_raw_reads=False,
    )

    assert spatial_table.exists()
    assert (tmp_path / "example_outputs" / "tables" / "qc_table.csv").exists()
