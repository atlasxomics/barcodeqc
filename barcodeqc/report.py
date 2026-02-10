from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import base64
import pandas as pd
import sys

from jinja2 import Template
from importlib.resources import files

import barcodeqc.utils as utils
import barcodeqc.paths as bq_paths
import barcodeqc.files as bq_files


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


def write_input_params(
    input_params: list[dict[str, str]],
    output_dir: Path,
    filename: str = "input_parameters.json",
) -> Path:
    import json

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / filename
    out_path.write_text(
        json.dumps(input_params, indent=2, ensure_ascii=True)
    )
    return out_path


def _image_data_uri(path: Path) -> str:
    with open(path, "rb") as image_file:
        encoded = base64.b64encode(image_file.read()).decode("utf-8")
    return f"data:image/png;base64,{encoded}"


def _figure_html(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".png", ".jpg", ".jpeg", ".gif", ".svg"}:
        data_uri = _image_data_uri(path)
        return f'<img src="{data_uri}" alt="{path.name}">'
    if suffix == ".html":
        return path.read_text()
    return f'<div class="note">Unsupported figure: {path.name}</div>'


def _barcode_set_from_params(
    input_params: Optional[list[dict[str, str]]],
) -> Optional[str]:
    if not input_params:
        return None
    for item in input_params:
        label = item.get("label")
        if label and label.startswith("_"):
            continue
        if label == "Barcode File":
            return item.get("value")
    return None


def _tissue_positions_provided(
    input_params: Optional[list[dict[str, str]]],
) -> Optional[bool]:
    if not input_params:
        return None
    for item in input_params:
        if item.get("label") == "_tissue_positions_provided":
            val = str(item.get("value", "")).strip().lower()
            if val in {"true", "1", "yes"}:
                return True
            if val in {"false", "0", "no"}:
                return False
    return None


def _top_n_by_label_from_barcode_set(
    barcode_set: Optional[str],
) -> Optional[dict[str, int]]:
    if not barcode_set:
        return None
    paths = bq_paths.BARCODE_PATHS.get(barcode_set)
    if not paths:
        return None
    try:
        bca = bq_files.open_barcode_file(paths["bca"])
        bcb = bq_files.open_barcode_file(paths["bcb"])
    except Exception:
        return None
    return {"L1": len(bca), "L2": len(bcb)}


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
    input_params: Optional[list[dict[str, str]]] = None,
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

    figure_list = list(figure_paths)
    by_stem: dict[str, Path] = {}
    for path in figure_list:
        stem = path.stem
        if stem not in by_stem:
            by_stem[stem] = path
            continue
        existing = by_stem[stem]
        if existing.suffix.lower() != ".html" and path.suffix.lower() == ".html":
            by_stem[stem] = path
    groups = _split_figures(by_stem.values())
    figures = {
        key: [_figure_html(_p) for _p in vals]
        for key, vals in groups.items()
    }
    bc_contam_figures = figures.get("bc_contam_l1", []) + figures.get(
        "bc_contam_l2", []
    )
    lane_qc_figures = figures.get("lane_qc_l1", []) + figures.get(
        "lane_qc_l2", []
    )

    display_params = None
    if input_params:
        display_params = [
            item for item in input_params
            if not str(item.get("label", "")).startswith("_")
        ]
    show_onoff_flag = _tissue_positions_provided(input_params)
    show_onoff = bool(onoff_table) if show_onoff_flag is None else show_onoff_flag

    bc_contam_labels: list[str] = []
    if figures["bc_contam_l1"]:
        bc_contam_labels.append("L1")
    if figures["bc_contam_l2"]:
        bc_contam_labels.append("L2")

    unexpected_barcodes: dict[str, list[dict[str, str]] | None] = {}
    unexpected_available = False
    unexpected_top_n = {"L1": 100, "L2": 100}
    if table_dir is not None and table_dir.exists():
        barcode_set = _barcode_set_from_params(input_params)
        top_n_by_label = _top_n_by_label_from_barcode_set(barcode_set)
        if top_n_by_label:
            unexpected_top_n.update(top_n_by_label)
        unexpected_barcodes, unexpected_available = (
            load_unexpected_barcodes_from_dir(
                table_dir,
                top_n_by_label=top_n_by_label,
            )
        )

    html_template = """
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>barcodeqc {{ sample_name }}</title>
  <link rel="stylesheet" href="{{ css_href }}">
  <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
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
            {% if show_onoff %}
            <li><a href="#onoff">On/Off Tissue</a></li>
            {% endif %}
          </ul>
        </nav>
      </aside>
      <main>

        <div class="intro">
          <p>
            This report summaries the main DBiT-seq failure modes upstream of
            alignment; test results are captured in the Summary table (below).
            Please review the table and supporting figures linked in the Table
            of Contents (left).
          </p>
        </div>

        {% if input_params %}
        <div class="input-metrics">
          <div class="input-metrics-grid">
            {% for item in input_params %}
            <div class="metric-card">
              <div class="metric-value">{{ item.value }}</div>
              <div class="metric-label">{{ item.label }}</div>
            </div>
            {% endfor %}
          </div>
        </div>
        {% else %}
        <div class="note">Input parameters not available.</div>
        {% endif %}

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
              {% set status_key = row.status | lower | replace(' ', '-') %}
              <td class="status {{ status_key }}">{{ row.status }}</td>
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
          {% if unexpected_available %}
          {% if bc_contam_labels %}
          <div class="barcode-table-carousel">
            {% for label in bc_contam_labels %}
            {% set rows = unexpected_barcodes.get(label) %}
            <div class="carousel-table{% if not loop.first %} hidden{% endif %}" data-carousel="bc-contam">
              <div class="table-title">{{ label }} unexpected barcodes in top {{ unexpected_top_n.get(label, 100) }} barcodes</div>
              {% if rows is none %}
              <div class="note">Unexpected barcode table not available.</div>
              {% else %}
              <div class="scroll-table">
              <table class="barcode-unexpected">
                <thead>
                  <tr>
                    <th>Barcode</th>
                    <th>Reads</th>
                    <th>Fraction</th>
                  </tr>
                </thead>
                <tbody>
                  {% if rows %}
                  {% for row in rows %}
                  <tr>
                    <td>{{ row.sequence }}</td>
                    <td>{{ row.count }}</td>
                    <td>{{ row.frac_count }}</td>
                  </tr>
                  {% endfor %}
                  {% else %}
                  <tr>
                    <td colspan="3">No unexpected barcodes found in the top {{ unexpected_top_n.get(label, 100) }}.</td>
                  </tr>
                  {% endif %}
                </tbody>
              </table>
              </div>
              {% endif %}
            </div>
            {% endfor %}
          </div>
          {% else %}
          <div class="note">Unexpected barcode table not available.</div>
          {% endif %}
          {% else %}
          <div class="note">Unexpected barcode table not available.</div>
          {% endif %}

          {% if bc_contam_figures %}
          <div class="plot-row">
            <div class="carousel-stage">
              {% for fig in bc_contam_figures %}
              <div class="carousel-figure{% if not loop.first %} hidden{% endif %}" data-carousel="bc-contam">
                {{ fig | safe }}
              </div>
              {% endfor %}
              {% if bc_contam_figures|length > 1 %}
              <button class="carousel-arrow prev" data-carousel="bc-contam" data-dir="-1">◀</button>
              <button class="carousel-arrow next" data-carousel="bc-contam" data-dir="1">▶</button>
              {% endif %}
            </div>
          </div>
          {% endif %}
        </div>

        <div id="lane-qc" class="row">
          <h2>Barcode Lane QC</h2>
          {% if lane_qc_figures %}
          <div class="plot-row">
            <div class="carousel-stage">
              {% for fig in lane_qc_figures %}
              <div class="carousel-figure{% if not loop.first %} hidden{% endif %}" data-carousel="lane-qc">
                {{ fig | safe }}
              </div>
              {% endfor %}
              {% if lane_qc_figures|length > 1 %}
              <button class="carousel-arrow prev" data-carousel="lane-qc" data-dir="-1">◀</button>
              <button class="carousel-arrow next" data-carousel="lane-qc" data-dir="1">▶</button>
              {% endif %}
            </div>
          </div>
          {% endif %}
        </div>

        {% if show_onoff %}
        <div id="onoff" class="row">
          <h2>On/Off Tissue Distribution</h2>

          <div class="onoff-grid">
            <div class="onoff-table">
              {% if onoff_table %}
              <table>
                {% for row in onoff_table %}
                <tr>
                  <td>{{ row.metric }}</td>
                  <td>{{ row.value }}</td>
                </tr>
                {% endfor %}
              </table>
              {% else %}
              <div class="note">On/off tissue metrics not available.</div>
              {% endif %}
            </div>

            <div class="onoff-plot">
              {% if figures.onoff %}
              <div class="plotly-figure">
                {{ figures.onoff[0] | safe }}
              </div>
              {% endif %}
            </div>

          </div>
        </div>
        {% endif %}

        {% if note_html %}
        <div class="note">{{ note_html | safe }}</div>
        {% endif %}
      </main>
    </div>
  </div>
  <script>
    function updateCarousel(name, dir) {
      const items = document.querySelectorAll(`.carousel-figure[data-carousel="${name}"]`);
      if (!items.length) return;
      const host = document.documentElement;
      const idxKey = `data-${name}-idx`;
      const cur = Number(host.getAttribute(idxKey) || 0);
      const next = (cur + dir + items.length) % items.length;
      items.forEach((item, idx) => {
        item.classList.toggle('hidden', idx !== next);
      });
      host.setAttribute(idxKey, String(next));
      const visible = items[next].querySelectorAll('.plotly-graph-div');
      if (window.Plotly && visible.length) {
        setTimeout(() => {
          visible.forEach((el) => Plotly.Plots.resize(el));
        }, 0);
      }
      const tables = document.querySelectorAll(`.carousel-table[data-carousel="${name}"]`);
      if (tables.length) {
        tables.forEach((table, idx) => {
          table.classList.toggle('hidden', idx !== next);
        });
      }
    }

    document.addEventListener('click', (e) => {
      const btn = e.target.closest('button[data-carousel]');
      if (btn) {
        const name = btn.getAttribute('data-carousel');
        const dir = Number(btn.getAttribute('data-dir') || 0);
        updateCarousel(name, dir);
      }
    });
    window.addEventListener('load', () => {
      if (!window.Plotly) return;
      document.querySelectorAll('.plotly-graph-div').forEach((el) => {
        Plotly.Plots.resize(el);
      });
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
        input_params=display_params,
        linker_metrics=linker_metrics,
        logo_uri=logo_uri,
        unexpected_barcodes=unexpected_barcodes,
        unexpected_available=unexpected_available,
        bc_contam_labels=bc_contam_labels,
        bc_contam_figures=bc_contam_figures,
        lane_qc_figures=lane_qc_figures,
        unexpected_top_n=unexpected_top_n,
        show_onoff=show_onoff,
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
            "Percent Pass Filtering": (
                f"{(adapter_reads / total_reads):.1%}"
                if total_reads > 0
                else "NA"
            ),
            "Number of Unique Barcodes": len(count_df),
            "Number Barcodes with 90% of reads": num_to_ninety,
            "Percent reads in expected barcodes": f"{total_read_from_expected:.1%}",
        }

    return metrics or None


def load_input_params_from_dir(
    sample_dir: Path,
    filename: str = "input_parameters.json",
) -> Optional[list[dict[str, str]]]:
    tables_dir = sample_dir / "tables"
    json_path = tables_dir / filename
    if not json_path.exists():
        json_path = sample_dir / filename
    if not json_path.exists():
        return None
    import json

    try:
        data = json.loads(json_path.read_text())
    except json.JSONDecodeError:
        return None
    if not isinstance(data, list):
        return None
    cleaned: list[dict[str, str]] = []
    for item in data:
        if not isinstance(item, dict):
            continue
        label = str(item.get("label", "")).strip()
        value = str(item.get("value", "")).strip()
        if not label and not value:
            continue
        cleaned.append({"label": label, "value": value})
    return cleaned or None


def load_unexpected_barcodes_from_dir(
    tables_dir: Path,
    top_n: int = 100,
    top_n_by_label: Optional[dict[str, int]] = None,
) -> tuple[dict[str, list[dict[str, str]] | None], bool]:
    rows: dict[str, list[dict[str, str]] | None] = {"L1": None, "L2": None}
    available = False
    for label in ("L1", "L2"):
        count_path = tables_dir / f"{label}_counts.csv"
        if not count_path.exists():
            continue

        count_df = pd.read_csv(count_path)
        if "expectMer" not in count_df.columns:
            continue

        available = True
        rows[label] = []
        sort_col = None
        if "count" in count_df.columns:
            sort_col = "count"
        elif "frac_count" in count_df.columns:
            sort_col = "frac_count"
        if sort_col:
            count_df = count_df.sort_values(sort_col, ascending=False)

        label_top_n = top_n_by_label.get(label, top_n) if top_n_by_label else top_n
        top = count_df.head(label_top_n).copy()
        expected_mask = top["expectMer"].astype(str).str.lower().isin(
            ["true", "1", "t", "yes"]
        )
        unexpected = top.loc[~expected_mask]
        if unexpected.empty:
            continue

        if "sequence" in unexpected.columns:
            seq_col = "sequence"
        else:
            fallback_cols = [
                c for c in unexpected.columns if not c.startswith("Unnamed")
            ]
            seq_col = fallback_cols[0] if fallback_cols else unexpected.columns[0]

        for _, row in unexpected.iterrows():
            count_val = row["count"] if "count" in unexpected.columns else None
            frac_val = (
                row["frac_count"]
                if "frac_count" in unexpected.columns
                else None
            )
            if pd.notna(count_val):
                try:
                    count_str = f"{int(count_val):,}"
                except (TypeError, ValueError):
                    count_str = str(count_val)
            else:
                count_str = "NA"
            if pd.notna(frac_val):
                try:
                    frac_str = f"{float(frac_val):.2%}"
                except (TypeError, ValueError):
                    frac_str = str(frac_val)
            else:
                frac_str = "NA"

            rows[label].append(
                {
                    "sequence": str(row[seq_col]),
                    "count": count_str,
                    "frac_count": frac_str,
                }
            )

    return rows, available
