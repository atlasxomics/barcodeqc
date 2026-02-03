from __future__ import annotations

import logging
import re

from pathlib import Path
from typing import List

logger = logging.getLogger(__name__)


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
