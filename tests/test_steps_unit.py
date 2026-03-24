from __future__ import annotations

from pathlib import Path

import pandas as pd

from barcodeqc.steps import (
    barcode_check_status,
    build_count_table,
    compute_hi_lo_qc,
    compute_onoff_metrics,
    lane_status,
    linker_conservation_status,
    make_spatial_table,
)


def _write_wc_file(path: Path, rows: list[tuple[str, str]]) -> None:
    path.write_text(
        "\n".join(f"x {barcode} {read_name}" for barcode, read_name in rows) + "\n",
        encoding="utf-8",
    )


def test_build_count_table_computes_expected_barcodes(tmp_path: Path) -> None:
    wc_path = tmp_path / "wc.txt"
    _write_wc_file(
        wc_path,
        [
            ("AAAACCCC", "read1"),
            ("AAAACCCC", "read2"),
            ("AAAACCCC", "read3"),
            ("TTTTGGGG", "read4"),
            ("TTTTGGGG", "read5"),
            ("CCCCAAAA", "read6"),
        ],
    )
    whitelist = pd.DataFrame(
        {
            "sequence": ["AAAACCCC", "TTTTGGGG"],
            "row": [1, 2],
            "col": [10, 20],
        }
    )

    count_table, unique_counts, expected_bcs, num_to_ninety, _ = (
        build_count_table(wc_path, whitelist, "row")
    )

    assert unique_counts["AAAACCCC"] == 3
    assert expected_bcs == {"AAAACCCC", "TTTTGGGG"}
    assert num_to_ninety == 2
    assert count_table.loc["AAAACCCC", "channel"] == 1
    assert bool(count_table.loc["CCCCAAAA", "expectMer"]) is False


def test_compute_hi_lo_qc_and_lane_status() -> None:
    count_table = pd.DataFrame(
        {
            "channel": [1, 2, 3, 4],
            "frac_count": [0.1, 0.1, 0.1, 0.7],
            "expectMer": [True, True, True, True],
        },
        index=["A", "B", "C", "D"],
    )
    count_table["sequence"] = count_table.index

    bc_table, total_hi_warn, total_lo_warn, total_mers = compute_hi_lo_qc(
        count_table
    )

    assert total_hi_warn == 1
    assert total_lo_warn == 3
    assert total_mers == 4
    assert lane_status(bc_table, "hiWarn") == "CAUTION"
    assert lane_status(bc_table.assign(hiWarn=False), "hiWarn") == "PASS"


def test_compute_onoff_metrics_returns_expected_values() -> None:
    spatial_table = pd.DataFrame(
        {
            "count": [10, 20, 5],
            "on_off": [1, 1, 0],
        }
    )

    metrics = compute_onoff_metrics(spatial_table).set_index("metric")["value"]

    assert metrics["total_pix"] == 3
    assert metrics["counts_on"] == 30
    assert metrics["counts_off"] == 5
    assert metrics["ratio_off_on"] == 5 / 30


def test_linker_conservation_status_uses_threshold() -> None:
    assert linker_conservation_status(100, 70) == ("PASS", 0.7)
    assert linker_conservation_status(100, 69) == ("CAUTION", 0.69)
    assert linker_conservation_status(0, 0) == ("CAUTION", 0.0)


def test_barcode_check_status_looks_at_top_n() -> None:
    count_table = pd.DataFrame(
        {
            "frac_count": [0.7, 0.2, 0.1],
            "expectMer": [True, False, True],
        },
        index=["A", "B", "C"],
    )

    assert barcode_check_status(count_table, top_n=1) == ("PASS", 0)
    assert barcode_check_status(count_table, top_n=2) == ("CAUTION", 1)


def test_make_spatial_table_merges_wildcards_and_positions(
    tmp_path: Path,
) -> None:
    wc1 = tmp_path / "wc1.txt"
    wc2 = tmp_path / "wc2.txt"
    positions = tmp_path / "positions.csv"

    _write_wc_file(
        wc1,
        [
            ("AAAACCCC", "read1"),
            ("AAAACCCC", "read2"),
            ("GGGGTTTT", "read3"),
            ("AAAACCCC", "read4"),
            ("AAAACCCC", "read5"),
            ("GGGGTTTT", "read6"),
        ],
    )
    _write_wc_file(
        wc2,
        [
            ("TTTTGGGG", "read1"),
            ("TTTTGGGG", "read2"),
            ("CCCCAAAA", "read3"),
            ("TTTTGGGG", "read4"),
            ("TTTTGGGG", "read5"),
            ("CCCCAAAA", "read6"),
        ],
    )
    positions.write_text(
        "\n".join(
            [
                "TTTTGGGGAAAACCCC,1,0,1",
                "CCCCAAAAGGGGTTTT,0,1,2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    spatial_table = make_spatial_table(wc1, wc2, positions)

    assert set(spatial_table["16mer"]) == {
        "TTTTGGGGAAAACCCC",
        "CCCCAAAAGGGGTTTT",
    }
    assert spatial_table["count"].sum() == 6
