from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def example_output_dir(repo_root: Path) -> Path:
    path = repo_root / "tests" / "data" / "example_outputs"
    if not path.exists():
        pytest.skip("tests/data example output fixture not present")
    return path


@pytest.fixture(scope="session")
def example_fastq_path(repo_root: Path) -> Path:
    path = (
        repo_root
        / "tests"
        / "data"
        / "example_fastq"
        / "barcodeqc_example_R2_001.fastq.gz"
    )
    if not path.exists():
        pytest.skip("tests/data example FASTQ fixture not present")
    return path


@pytest.fixture(scope="session")
def example_tissue_positions_path(repo_root: Path) -> Path:
    path = (
        repo_root
        / "tests"
        / "data"
        / "example_fastq"
        / "tissue_positions_list.csv"
    )
    if not path.exists():
        pytest.skip("tests/data example tissue positions fixture not present")
    return path
