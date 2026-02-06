from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import base64
import pandas as pd
import sys

from jinja2 import Template
from importlib.resources import files

import barcodeqc.utils as utils


def _load_static_logo() -> str | None:
    logo_path = files("barcodeqc") / "data" / "static" / "logo.png"
    if not logo_path.is_file():
        return None
    return _image_data_uri(logo_path)


def write_summary_table(
    summary_table: pd.DataFrame,
    output_dir: Path,
    filename: str = "qc_table.csv",
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / filename
    summary_table.to_csv(out_path, index=False)
    return out_path


def _image_data_uri(path: Path) -> str:
    with open(path, "rb") as image_file:
        encoded = base64.b64encode(image_file.read()).decode("utf-8")
    return f"data:image/png;base64,{encoded}"


def _split_figures(figure_paths: Iterable[Path]) -> dict[str, list[Path]]:
    groups: dict[str, list[Path]] = {
        "bc_contam_l1": [],
        "bc_contam_l2": [],
        "lane_qc_l1": [],
        "lane_qc_l2": [],
        "onoff": [],
        "other": [],
    }
    for p in figure_paths:
        name = p.name.lower()
        if "pareto" in name and "l1" in name:
            groups["bc_contam_l1"].append(p)
        elif "pareto" in name and "l2" in name:
            groups["bc_contam_l2"].append(p)
        elif "barplot" in name and "l1" in name:
            groups["lane_qc_l1"].append(p)
        elif "barplot" in name and "l2" in name:
            groups["lane_qc_l2"].append(p)
        elif "dense_on_off" in name:
            groups["onoff"].append(p)
        else:
            groups["other"].append(p)
    return groups


def print_summary_table(summary_table: pd.DataFrame) -> None:
    if summary_table.empty:
        return

    headers = ["METRIC", "STATUS"]
    rows = summary_table[["metric", "status"]].astype(str).values.tolist()
    widths = [
        max(len(headers[0]), *(len(r[0]) for r in rows)),
        max(len(headers[1]), *(len(r[1]) for r in rows)),
    ]

    def fmt_row(left: str, right: str) -> str:
        return f"| {left.ljust(widths[0])} | {right.ljust(widths[1])} |"

    border = f"+-{'-' * widths[0]}-+-{'-' * widths[1]}-+"

    print(border, file=sys.stdout)
    print(fmt_row(headers[0], headers[1]), file=sys.stdout)
    print(border, file=sys.stdout)
    for left, right in rows:
        print(fmt_row(left, right), file=sys.stdout)
    print(border, file=sys.stdout)


def generate_report(
    figure_paths: Iterable[Path],
    output_dir: Path,
    sample_name: str,
    note_html: str = "",
    summary_table: Optional[pd.DataFrame] = None,
    linker_metrics: Optional[dict[str, dict[str, str | int | float]]] = None,
    onoff_table: Optional[pd.DataFrame] = None,
    file_tag: str = "bcQC",
    table_dir: Optional[Path] = None,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)

    logo_uri = _load_static_logo()

    outputs: dict[str, Path] = {}
    if summary_table is not None:
        target_dir = table_dir or output_dir
        outputs["summary_csv"] = write_summary_table(
            summary_table,
            target_dir,
        )

    groups = _split_figures(figure_paths)
    figures = {
        key: [(_p.name, _image_data_uri(_p)) for _p in vals]
        for key, vals in groups.items()
    }

    html_template = """
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>barcodeqc {{ sample_name }}</title>
  <link rel="stylesheet" href="{{ css_href }}">
</head>
<body>
  <div class="container">

    <div class="topbar">
      <div class="topbar-inner">

        <div class="topbar-logo">
          {% if logo_uri %}
            <img src="{{ logo_uri }}" alt="logo" class="logo">
          {% endif %}
        </div>

        <div class="topbar-title">
          <h1>Barcode QC Report</h1>
        </div>

        <div class="topbar-pill">
          <div class="sample-pill">{{ sample_name }}</div>
        </div>

      </div>
    </div>

    <div class="layout">

      <aside class="sidebar">
        <nav class="nav">
          <ul class="nav-list">
            <li><a href="#summary">Summary</a></li>
            <li><a href="#linker-filtering">Linker Filtering</a></li>
            <li><a href="#barcode-check">Barcode Check</a></li>
            <li><a href="#lane-qc">Lane QC</a></li>
            <li><a href="#onoff">On/Off Tissue</a></li>
          </ul>
        </nav>
      </aside>
      <main>

        <div class="intro">
          <p>
            `barcodeqc` is a lightweight command-line tool for the rapid evaluation
            of AtlasXomics epigenomic DBiT-seq experiments. The tool takes fastq
            files from read2 an input and analyzes the distribution of barcode
            sequences to screen for the top preliminary failure modes (low linker
            conservation, mismatched barcodes, missing/overabundent barcodes).
          </p>
        </div>

        <div id="summary" class="summary-panel row">
          <div class="summary-content">
          <h2>Summary</h2>
          {% if summary_table %}
          <table class="summary-table">
            <colgroup>
              <col style="width: 15%">
              <col style="width: 15%">
              <col style="width: 70%">
            </colgroup>
            {% for row in summary_table %}
            <tr>
              <td>{{ row.metric }}</td>
              <td>{{ row.status }}</td>
              <td>{{ row.description }}</td>
            </tr>
            {% endfor %}
          </table>
          {% else %}
          <div class="note">Summary table not available.</div>
          {% endif %}
          </div>
        </div>

        <div id="linker-filtering" class="row">
          <h2>Linker Filtering</h2>
          {% if linker_metrics %}
          <div class="grid-2">
            {% for label, metrics in linker_metrics.items() %}
            <table>
              {% for key, value in metrics.items() %}
              <tr><td>{{ label }} {{ key }}</td><td>{{ value }}</td></tr>
              {% endfor %}
            </table>
            {% endfor %}
          </div>
          {% else %}
          <div class="note">Linker metrics not available.</div>
          {% endif %}
        </div>

        <div id="barcode-check" class="row">
          <h2>Barcode Check</h2>
          {% if figures.bc_contam_l1 %}
          <div class="plot-row">
            <div class="carousel-stage">
              <img id="bc-contam-img" src="{{ figures.bc_contam_l1[0][1] }}" alt="{{ figures.bc_contam_l1[0][0] }}">
              {% if figures.bc_contam_l2 %}
              <button class="carousel-arrow prev" data-carousel="bc-contam" data-dir="-1">◀</button>
              <button class="carousel-arrow next" data-carousel="bc-contam" data-dir="1">▶</button>
              {% endif %}
            </div>
          </div>
          {% endif %}
        </div>

        <div id="lane-qc" class="row">
          <h2>Barcode Lane QC</h2>
          {% if figures.lane_qc_l1 %}
          <div class="plot-row">
            <div class="carousel-stage">
              <img id="lane-qc-img" src="{{ figures.lane_qc_l1[0][1] }}" alt="{{ figures.lane_qc_l1[0][0] }}">
              {% if figures.lane_qc_l2 %}
              <button class="carousel-arrow prev" data-carousel="lane-qc" data-dir="-1">◀</button>
              <button class="carousel-arrow next" data-carousel="lane-qc" data-dir="1">▶</button>
              {% endif %}
            </div>
          </div>
          {% endif %}
        </div>

        <div id="onoff" class="row">
          <h2>On/Off Tissue Distribution</h2>
          <div class="narrow">
            {% if onoff_table %}
            <table>
              {% for row in onoff_table %}
              <tr><td>{{ row.metric }}</td><td>{{ row.value }}</td></tr>
              {% endfor %}
            </table>
            {% else %}
            <div class="note">On/off tissue metrics not available.</div>
            {% endif %}
            <div class="onoff-plot">
              {% if figures.onoff %}
              <img src="{{ figures.onoff[0][1] }}" alt="{{ figures.onoff[0][0] }}">
              {% endif %}
            </div>
          </div>
        </div>

        {% if note_html %}
        <div class="note">{{ note_html | safe }}</div>
        {% endif %}
      </main>
    </div>
  </div>
  <script>
    const carousels = {
      "bc-contam": [
        {% if figures.bc_contam_l1 %}"{{ figures.bc_contam_l1[0][1] }}",{% endif %}
        {% if figures.bc_contam_l2 %}"{{ figures.bc_contam_l2[0][1] }}"{% endif %}
      ],
      "lane-qc": [
        {% if figures.lane_qc_l1 %}"{{ figures.lane_qc_l1[0][1] }}",{% endif %}
        {% if figures.lane_qc_l2 %}"{{ figures.lane_qc_l2[0][1] }}"{% endif %}
      ]
    };

    function updateCarousel(name, dir) {
      const imgs = carousels[name] || [];
      if (!imgs.length) return;
      const imgEl = document.getElementById(`${name}-img`);
      const idxKey = `data-${name}-idx`;
      const cur = Number(imgEl.getAttribute(idxKey) || 0);
      const next = (cur + dir + imgs.length) % imgs.length;
      imgEl.src = imgs[next];
      imgEl.setAttribute(idxKey, String(next));
    }

    document.addEventListener('click', (e) => {
      const btn = e.target.closest('button[data-carousel]');
      if (btn) {
        const name = btn.getAttribute('data-carousel');
        const dir = Number(btn.getAttribute('data-dir') || 0);
        updateCarousel(name, dir);
      }
    });
  </script>
</body>
</html>
    """

    template = Template(html_template)
    html_content = template.render(
        sample_name=sample_name,
        css_href=files("barcodeqc") / "data" / "static" / "report.css",
        note_html=note_html,
        figures=figures,
        summary_table=(
            summary_table.to_dict(orient="records")
            if summary_table is not None else None
        ),
        onoff_table=(
            onoff_table.to_dict(orient="records")
            if onoff_table is not None else None
        ),
        linker_metrics=linker_metrics,
        logo_uri=logo_uri,
    )

    html_path = output_dir / f"{sample_name}_{file_tag}_report.html"
    html_path.write_text(html_content)
    outputs["html"] = html_path

    return outputs


def load_linker_metrics_from_dir(
    sample_dir: Path,
) -> Optional[dict[str, dict[str, str | int | float]]]:
    metrics: dict[str, dict[str, str | int | float]] = {}
    logs_dir = sample_dir / "logs"
    tables_dir = sample_dir / "tables"
    for label in ("L1", "L2"):
        log_path = logs_dir / f"cutadapt_{label}.log"
        count_path = tables_dir / f"{label}_counts.csv"
        if not log_path.exists():
            log_path = sample_dir / f"cutadapt_{label}.log"
        if not count_path.exists():
            count_path = sample_dir / f"{label}_counts.csv"
        if not log_path.exists() or not count_path.exists():
            continue

        total_reads_str, adapter_reads_str = utils.parse_read_log(log_path)
        total_reads = int(total_reads_str)
        adapter_reads = int(adapter_reads_str)

        count_df = pd.read_csv(count_path)
        if "frac_count" not in count_df.columns:
            continue

        if "cumulative_sum" in count_df.columns:
            num_to_ninety = int((count_df["cumulative_sum"] <= 0.9).sum())
        else:
            sorted_frac = count_df["frac_count"].sort_values(ascending=False)
            num_to_ninety = int((sorted_frac.cumsum() <= 0.9).sum())

        frac_sorted = count_df["frac_count"].sort_values(ascending=False)
        pct_for_50 = f"{frac_sorted.iloc[:50].sum():.1%}"
        pct_for_96 = f"{frac_sorted.iloc[:96].sum():.1%}"

        if "expectMer" in count_df.columns:
            expected_mask = count_df["expectMer"].astype(str).str.lower().isin(
                ["true", "1"]
            )
            total_read_from_expected = count_df.loc[
                expected_mask, "frac_count"
            ].sum()
        else:
            total_read_from_expected = 0.0

        metrics[label] = {
            "Total Reads": total_reads,
            "Total with Linker": adapter_reads,
            "Number of Unique Barcodes": len(count_df),
            "Number Barcodes with 90% of reads": num_to_ninety,
            "Percent reads in expected barcodes": f"{total_read_from_expected:.1%}",
            "Percent reads in top 50 8mers": pct_for_50,
            "Percent reads in top 96 8mers": pct_for_96,
        }

    return metrics or None
