from __future__ import annotations

from pathlib import Path

import pytest

import barcodeqc.paths as paths
from barcodeqc.files import (
    BarcodeFileError,
    WildcardFileError,
    load_wc_file,
    open_barcode_file,
    open_positions_file,
)


def test_open_barcode_file_accepts_builtin_fixture() -> None:
    barcodes = open_barcode_file(paths.BARCODE_PATHS["bc50"]["bca"])

    assert list(barcodes.columns) == ["sequence", "row", "col"]
    assert not barcodes.empty


def test_open_barcode_file_rejects_bad_sequence(tmp_path: Path) -> None:
    path = tmp_path / "bad_barcodes.csv"
    path.write_text(
        "sequence,row,col\nBADSEQ,1,2\n",
        encoding="utf-8",
    )

    with pytest.raises(BarcodeFileError, match="Invalid barcodes detected"):
        open_barcode_file(path)


def test_open_positions_file_strips_trailing_dash_one(tmp_path: Path) -> None:
    path = tmp_path / "positions.csv"
    path.write_text(
        "AAAACCCCGGGGTTTT-1,1,0,1\nTTTTGGGGCCCCAAAA,0,1,2\n",
        encoding="utf-8",
    )

    positions = open_positions_file(path)

    assert positions["barcodes"].tolist() == [
        "AAAACCCCGGGGTTTT",
        "TTTTGGGGCCCCAAAA",
    ]
    assert positions["on_off"].tolist() == [1, 0]


def test_open_positions_file_rejects_invalid_on_off(tmp_path: Path) -> None:
    path = tmp_path / "positions.csv"
    path.write_text(
        "AAAACCCCGGGGTTTT,2,0,1\n",
        encoding="utf-8",
    )

    with pytest.raises(BarcodeFileError, match="on_off column must contain only 0 or 1"):
        open_positions_file(path)


def test_load_wc_file_accepts_valid_fixture(tmp_path: Path) -> None:
    path = tmp_path / "wc.txt"
    path.write_text(
        "\n".join(
            [
                "AAAACCCC read1",
                "TTTTGGGG read2",
                "AAAACCCC read3",
                "TTTTGGGG read4",
                "AAAACCCC read5",
                "TTTTGGGG read6",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    wc = load_wc_file(path)

    assert wc["8mer"].tolist()[:2] == ["AAAACCCC", "TTTTGGGG"]
    assert wc["readName"].tolist()[:2] == ["read1", "read2"]


def test_load_wc_file_accepts_variable_whitespace(tmp_path: Path) -> None:
    path = tmp_path / "wc.txt"
    path.write_text(
        "\n".join(
            [
                "AAAACCCC   read1",
                "TTTTGGGG    read2",
                "AAAACCCC read3",
                "TTTTGGGG read4",
                "AAAACCCC   read5",
                "TTTTGGGG read6",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    wc = load_wc_file(path)

    assert wc["8mer"].tolist()[:2] == ["AAAACCCC", "TTTTGGGG"]
    assert wc["readName"].tolist()[:2] == ["read1", "read2"]


def test_load_wc_file_preserves_empty_barcode_captures(tmp_path: Path) -> None:
    path = tmp_path / "wc.txt"
    path.write_text(
        "\n".join(
            [
                "  lh00134:1:1:1 1:N:0:CGTACTAG",
                "AGGTACTC lh00134:1:1:2 1:N:0:CGTACTAG",
                " lh00134:1:1:3 1:N:0:CGTACTAG",
                "TTTTGGGG lh00134:1:1:4 1:N:0:CGTACTAG",
                "AAAACCCC lh00134:1:1:5 1:N:0:CGTACTAG",
                " lh00134:1:1:6 1:N:0:CGTACTAG",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    wc = load_wc_file(path)

    assert wc["8mer"].tolist()[:4] == ["", "AGGTACTC", "", "TTTTGGGG"]
    assert wc["readName"].tolist()[:2] == [
        "lh00134:1:1:1 1:N:0:CGTACTAG",
        "lh00134:1:1:2 1:N:0:CGTACTAG",
    ]


def test_load_wc_file_rejects_short_input(tmp_path: Path) -> None:
    path = tmp_path / "wc.txt"
    path.write_text(
        "\n".join(
            [
                "AAAACCCC read1",
                "TTTTGGGG read2",
                "AAAACCCC read3",
                "TTTTGGGG read4",
                "AAAACCCC read5",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(WildcardFileError, match="Fewer than 6 reads found"):
        load_wc_file(path)
