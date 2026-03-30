from __future__ import annotations

import logging
import re

from pathlib import Path
from shutil import which
from typing import List
import gzip

logger = logging.getLogger(__name__)


class ExternalDependencyError(RuntimeError):
    pass


class InputFastqError(ValueError):
    """Raised when the supplied FASTQ is incompatible with barcodeqc."""


def contains_acgt_word(input_list: List[str]) -> List[int]:
    '''Function to check for 8-character word made up of A, C, G, T in a list
    and return indices.'''
    pattern = re.compile(r'^[ACGTN]{8,}$')
    return [
        index
        for index, item in enumerate(input_list)
        if isinstance(item, str) and pattern.search(item)
    ]


def infer_fastq_read_number(path: Path) -> int | None:
    name = path.name
    normalized = name
    for suffix in (".gz", ".gzip", ".fastq", ".fq"):
        if normalized.lower().endswith(suffix):
            normalized = normalized[:-len(suffix)]

    explicit = re.search(r"(?:^|[._-])R([12])(?:$|[._-])", normalized, re.I)
    if explicit:
        return int(explicit.group(1))

    terminal = re.search(r"(?:^|[._-])([12])$", normalized)
    if terminal:
        return int(terminal.group(1))

    return None


def validate_r2_fastq_path(path: Path) -> None:
    inferred_read = infer_fastq_read_number(path)
    if inferred_read != 1:
        return

    logger.warning(
        f"Input FASTQ appears to be Read 1 based on its filename: {path.name}. "
        "barcodeqc requires the Read 2 FASTQ because AtlasXomics barcodes and "
        "ligation linkers are expected in Read 2. Please verify the input file "
        "and rerun with the Read 2 FASTQ, typically a file containing 'R2' in "
        "its name."
    )


def parse_read_log(log_path: str) -> tuple[str, str]:
    text = Path(log_path).read_text()

    tot = re.search(r"Total reads processed:\s*([\d,]+)", text)
    adapt = re.search(r"Reads with adapters:\s*([\d,]+)", text)

    if not tot or not adapt:
        raise ValueError("Expected read counts not found")

    total_reads = tot.group(1).replace(",", "")
    adapter_reads = adapt.group(1).replace(",", "")
    return total_reads, adapter_reads


def validate_linker_detection(
    r2_path: Path,
    log_linker1: Path,
    log_linker2: Path,
) -> None:
    _, adapter_reads_l1 = parse_read_log(log_linker1)
    _, adapter_reads_l2 = parse_read_log(log_linker2)

    if int(adapter_reads_l1) > 0 or int(adapter_reads_l2) > 0:
        return

    raise InputFastqError(
        f"No AtlasXomics linker sequences were detected in {r2_path.name}. "
        "Please confirm that you supplied the Read 2 FASTQ, not Read 1. "
        "barcodeqc expects Read 2 in the format "
        "'linker1 | barcodeA | linker2 | barcodeB | genomic sequence'."
    )


def require_executable(name: str) -> str:
    exe = which(name)
    if not exe:
        raise ExternalDependencyError(
            f"{name} not found on PATH. Please install {name} and try again."
        )
    return exe


def count_fastq_reads(path: Path) -> int:
    opener = gzip.open if path.suffix in {".gz", ".gzip"} else open
    line_count = 0
    with opener(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            line_count += chunk.count(b"\n")
    if line_count % 4 != 0:
        raise ValueError(
            f"Fastq file has {line_count} lines (not divisible by 4): {path}"
        )
    return line_count // 4
