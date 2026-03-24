from __future__ import annotations

from pathlib import Path

import pandas as pd

import barcodeqc.report as report


def test_generate_report_accepts_onoff_table_without_hidden_flag(
    tmp_path: Path,
) -> None:
    figure_path = tmp_path / "L1_pareto.html"
    figure_path.write_text("<div>plot</div>", encoding="utf-8")

    summary_table = pd.DataFrame(
        [{"metric": "Linker 1 Filter", "status": "PASS", "description": "ok"}]
    )
    onoff_table = pd.DataFrame(
        [{"metric": "ratio_off_on", "value": 0.1}]
    )

    outputs = report.generate_report(
        figure_paths=[figure_path],
        output_dir=tmp_path,
        sample_name="sample",
        summary_table=summary_table,
        onoff_table=onoff_table,
        input_params=[{"label": "Barcode File", "value": "bc50"}],
    )

    html = outputs["html"].read_text(encoding="utf-8")

    assert outputs["html"].exists()
    assert "On/Off Tissue" in html


def test_load_input_params_from_dir_reads_fixture(
    example_output_dir: Path,
) -> None:
    params = report.load_input_params_from_dir(example_output_dir)

    assert params is not None
    assert {"label": "Barcode File", "value": "bc220_20-MAY"} in params


def test_load_linker_metrics_from_dir_reads_fixture(
    example_output_dir: Path,
) -> None:
    metrics = report.load_linker_metrics_from_dir(example_output_dir)

    assert metrics is not None
    assert set(metrics) == {"L1", "L2"}
    assert metrics["L1"]["Total Reads"] > 0


def test_load_unexpected_barcodes_from_dir_reads_fixture(
    example_output_dir: Path,
) -> None:
    rows, available = report.load_unexpected_barcodes_from_dir(
        example_output_dir / "tables"
    )

    assert available is True
    assert rows["L1"] is not None
    assert rows["L2"] is not None
    assert any(row["sequence"] == "ACCGACAT" for row in rows["L1"])
