from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from barcodeqc.config import QCConfig
from barcodeqc.steps import run_cutadapt, run_subsample


class DummyPipe:
    def close(self) -> None:
        return None


class DummyPopen:
    def __init__(self, returncode: int = 0) -> None:
        self.stdout = DummyPipe()
        self._returncode = returncode

    def wait(self) -> int:
        return self._returncode


def test_run_subsample_prefers_pigz(monkeypatch, tmp_path: Path) -> None:
    calls: dict[str, list[str]] = {}

    def fake_popen(cmd, stdout=None):
        calls["seqtk"] = cmd
        return DummyPopen()

    def fake_run(cmd, stdin=None, stdout=None, check=None):
        calls["gzip"] = cmd
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("barcodeqc.steps.which", lambda name: "/usr/bin/pigz")
    monkeypatch.setattr("barcodeqc.steps.subprocess.Popen", fake_popen)
    monkeypatch.setattr("barcodeqc.steps.subprocess.run", fake_run)

    config = QCConfig(
        sample_name="sample",
        r2_path=tmp_path / "reads.fastq.gz",
        barcode_set="bc50",
        sample_reads=1000,
        random_seed=13,
        tissue_position_file=None,
        output_dir=tmp_path / "out",
    )

    ds_path = run_subsample(config, tmp_path)

    assert ds_path == tmp_path / "ds_1000.fastq.gz"
    assert ds_path.exists()
    assert calls["seqtk"][:4] == ["seqtk", "sample", "-s", "13"]
    assert calls["gzip"] == ["/usr/bin/pigz", "-1"]


def test_run_cutadapt_writes_logs(monkeypatch, tmp_path: Path) -> None:
    calls: list[list[str]] = []

    def fake_run(cmd, stdout=None, stderr=None, check=None):
        calls.append(cmd)
        stdout.write(
            "Total reads processed: 10\nReads with adapters: 8 (80.0%)\n"
        )
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("barcodeqc.steps.subprocess.run", fake_run)

    ds_path = tmp_path / "reads.fastq.gz"
    ds_path.write_text("", encoding="utf-8")

    wc1, log1, wc2, log2 = run_cutadapt(ds_path, tmp_path, cores=4)

    assert [cmd[0] for cmd in calls] == ["cutadapt", "cutadapt"]
    assert "--cores" in calls[0]
    assert "4" in calls[0]
    assert log1.exists()
    assert log2.exists()
    assert wc1.name == "cutadapt_wc_L1.txt"
    assert wc2.name == "cutadapt_wc_L2.txt"
