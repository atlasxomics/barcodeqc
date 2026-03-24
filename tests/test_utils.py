from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from barcodeqc.utils import contains_acgt_word, count_fastq_reads, parse_read_log


def test_contains_acgt_word_returns_matching_indices() -> None:
    assert contains_acgt_word(["x", "AAAACCCC", "nope", "TTTTGGGG"]) == [1, 3]


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
