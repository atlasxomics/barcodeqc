from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from barcodeqc.utils import (
    InputFastqError,
    contains_acgt_word,
    count_fastq_reads,
    infer_fastq_read_number,
    log_barcode_capture_quality,
    parse_read_log,
    validate_linker_detection,
)


def test_contains_acgt_word_returns_matching_indices() -> None:
    assert contains_acgt_word(["x", "AAAACCCC", "nope", "TTTTGGGG"]) == [1, 3]


def test_contains_acgt_word_ignores_non_string_values() -> None:
    assert contains_acgt_word(["x", float("nan"), "AAAACCCC"]) == [2]


def test_infer_fastq_read_number_detects_common_r1_r2_names() -> None:
    assert infer_fastq_read_number(Path("sample_R1_001.fastq.gz")) == 1
    assert infer_fastq_read_number(Path("sample_R2_001.fastq.gz")) == 2
    assert infer_fastq_read_number(Path("sample_1.fastq.gz")) == 1
    assert infer_fastq_read_number(Path("sample.fastq.gz")) is None


def test_parse_read_log_extracts_counts(tmp_path: Path) -> None:
    path = tmp_path / "cutadapt.log"
    path.write_text(
        "Total reads processed: 1,234\nReads with adapters: 567 (45.9%)\n",
        encoding="utf-8",
    )

    assert parse_read_log(path) == ("1234", "567")


def test_count_fastq_reads_handles_gzip(tmp_path: Path) -> None:
    path = tmp_path / "reads.fastq.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("@r1\nACGT\n+\n!!!!\n@r2\nTGCA\n+\n####\n")

    assert count_fastq_reads(path) == 2


def test_count_fastq_reads_rejects_invalid_line_count(tmp_path: Path) -> None:
    path = tmp_path / "bad.fastq"
    path.write_text("@r1\nACGT\n+\n", encoding="utf-8")

    with pytest.raises(ValueError, match="not divisible by 4"):
        count_fastq_reads(path)


def test_validate_linker_detection_rejects_missing_linkers(tmp_path: Path) -> None:
    log1 = tmp_path / "cutadapt_L1.log"
    log2 = tmp_path / "cutadapt_L2.log"
    log1.write_text(
        "Total reads processed: 100\nReads with adapters: 0 (0.0%)\n",
        encoding="utf-8",
    )
    log2.write_text(
        "Total reads processed: 100\nReads with adapters: 0 (0.0%)\n",
        encoding="utf-8",
    )

    with pytest.raises(InputFastqError, match="Read 2 FASTQ, not Read 1"):
        validate_linker_detection(
            Path("sample_R1_001.fastq.gz"),
            log1,
            log2,
        )


def test_log_barcode_capture_quality_logs_error_without_raising(
    caplog,
) -> None:
    log_barcode_capture_quality(
        "L1",
        adapter_reads=100,
        valid_barcode_reads=70,
        empty_barcode_reads=30,
    )

    assert "30.0%" in caplog.text
    assert "Continuing the run" in caplog.text
