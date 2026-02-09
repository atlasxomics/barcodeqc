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


def contains_acgt_word(input_list: List[str]) -> List[int]:
    '''Function to check for 8-character word made up of A, C, G, T in a list
    and return indices.'''
    pattern = re.compile(r'^[ACGTN]{8,}$')
    return [
        index for index, item in enumerate(input_list) if pattern.search(item)
    ]


def parse_read_log(log_path: str) -> tuple[str, str]:
    text = Path(log_path).read_text()

    tot = re.search(r"Total reads processed:\s*([\d,]+)", text)
    adapt = re.search(r"Reads with adapters:\s*([\d,]+)", text)

    if not tot or not adapt:
        raise ValueError("Expected read counts not found")

    total_reads = tot.group(1).replace(",", "")
    adapter_reads = adapt.group(1).replace(",", "")
    return total_reads, adapter_reads


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
