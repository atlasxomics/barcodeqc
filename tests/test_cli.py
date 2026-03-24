from __future__ import annotations

from pathlib import Path

import barcodeqc.cli as cli
from barcodeqc.utils import ExternalDependencyError


def test_parse_args_for_qc(tmp_path: Path) -> None:
    args = cli.parse_args(
        [
            "qc",
            "sample1",
            str(tmp_path / "reads.fastq.gz"),
            "bc50",
            "-r",
            "123",
            "-s",
            "7",
        ]
    )

    assert args.command == "qc"
    assert args.sample_name == "sample1"
    assert args.r2_path == tmp_path / "reads.fastq.gz"
    assert args.barcode_set == "bc50"
    assert args.sample_reads == 123
    assert args.random_seed == 7


def test_run_returns_usage_on_invalid_args() -> None:
    assert cli.run([]) == cli.EX_USAGE


def test_run_qc_invokes_pipeline(monkeypatch, tmp_path: Path) -> None:
    out_dir = tmp_path / "sample_outputs"
    spatial_table = out_dir / "tables" / "spatialTable.csv"

    def fake_qc(**kwargs):
        spatial_table.parent.mkdir(parents=True, exist_ok=True)
        spatial_table.write_text("count\n1\n", encoding="utf-8")
        return spatial_table

    monkeypatch.setattr(cli, "setup_logging", lambda **kwargs: None)
    monkeypatch.setattr(cli, "output_dir_from_sample_name", lambda _: out_dir)
    monkeypatch.setattr(cli, "qc", fake_qc)

    code = cli.run(
        ["qc", "sample", str(tmp_path / "reads.fastq.gz"), "bc50"]
    )

    assert code == cli.EX_OK
    assert spatial_table.exists()


def test_run_maps_external_dependency_error(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(cli, "setup_logging", lambda **kwargs: None)
    monkeypatch.setattr(
        cli,
        "output_dir_from_sample_name",
        lambda _: tmp_path / "sample_outputs",
    )
    monkeypatch.setattr(
        cli,
        "qc",
        lambda **kwargs: (_ for _ in ()).throw(
            ExternalDependencyError("seqtk missing")
        ),
    )

    code = cli.run(
        ["qc", "sample", str(tmp_path / "reads.fastq.gz"), "bc50"]
    )

    assert code == cli.EX_UNAVAILABLE


def test_report_command_generates_output_from_existing_run(
    monkeypatch,
    tmp_path: Path,
    example_output_dir: Path,
) -> None:
    out_dir = tmp_path / "report_outputs"

    monkeypatch.setattr(cli, "setup_logging", lambda **kwargs: None)
    monkeypatch.setattr(cli, "output_dir_from_sample_name", lambda _: out_dir)

    code = cli.run(["report", "example", str(example_output_dir)])

    assert code == cli.EX_OK
    assert (out_dir / "example_bcQC_report.html").exists()
