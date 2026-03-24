from __future__ import annotations

from pathlib import Path

import pytest

import barcodeqc.paths as paths
from barcodeqc.config import QCConfig, output_dir_from_sample_name


def test_output_dir_from_sample_name_appends_suffix(
    monkeypatch,
    tmp_path: Path,
) -> None:
    monkeypatch.chdir(tmp_path)

    assert output_dir_from_sample_name("sample") == tmp_path / "sample_outputs"
    assert output_dir_from_sample_name("sample_outputs") == (
        tmp_path / "sample_outputs"
    )


def test_qcconfig_from_args_uses_default_positions(
    monkeypatch,
    tmp_path: Path,
) -> None:
    monkeypatch.chdir(tmp_path)

    config = QCConfig.from_args(
        sample_name="sample",
        r2_path=tmp_path / "reads.fastq.gz",
        barcode_set="bc50",
        sample_reads=1000,
        random_seed=42,
        tissue_position_file=None,
    )

    assert config.output_dir == tmp_path / "sample_outputs"
    assert config.tissue_position_file == paths.BARCODE_PATHS["bc50"]["positions"]


def test_qcconfig_validate_checks_required_files(
    monkeypatch,
    tmp_path: Path,
) -> None:
    r2_path = tmp_path / "reads.fastq.gz"
    position_path = tmp_path / "positions.csv"
    r2_path.write_text("", encoding="utf-8")
    position_path.write_text("", encoding="utf-8")

    monkeypatch.setattr("barcodeqc.config.require_executable", lambda _: "seqtk")

    config = QCConfig(
        sample_name="sample",
        r2_path=r2_path,
        barcode_set="bc50",
        sample_reads=1000,
        random_seed=42,
        tissue_position_file=position_path,
        output_dir=tmp_path / "out",
    )

    config.validate()


def test_qcconfig_validate_raises_for_missing_fastq(tmp_path: Path) -> None:
    config = QCConfig(
        sample_name="sample",
        r2_path=tmp_path / "missing.fastq.gz",
        barcode_set="bc50",
        sample_reads=1000,
        random_seed=42,
        tissue_position_file=None,
        output_dir=tmp_path / "out",
    )

    with pytest.raises(FileNotFoundError, match="fastq file path does not exist"):
        config.validate()


def test_qcconfig_validate_raises_for_missing_positions(
    monkeypatch,
    tmp_path: Path,
) -> None:
    r2_path = tmp_path / "reads.fastq.gz"
    r2_path.write_text("", encoding="utf-8")

    config = QCConfig(
        sample_name="sample",
        r2_path=r2_path,
        barcode_set="bc50",
        sample_reads=1000,
        random_seed=42,
        tissue_position_file=tmp_path / "missing_positions.csv",
        output_dir=tmp_path / "out",
    )

    monkeypatch.setattr("barcodeqc.config.require_executable", lambda _: "seqtk")

    with pytest.raises(FileNotFoundError, match="Could not find tissue_postion file"):
        config.validate()
