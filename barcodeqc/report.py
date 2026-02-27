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


def _load_static_image(filename: str) -> str | None:
    image_path = files("barcodeqc") / "data" / "static" / filename
    if not image_path.is_file():
        return None
    return _image_data_uri(image_path)


def _load_static_text(filename: str) -> str | None:
    text_path = files("barcodeqc") / "data" / "static" / filename
    if not text_path.is_file():
        return None
    return text_path.read_text(encoding="utf-8")


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

    logo_uri = _load_static_image("logo.png")
    css_text = _load_static_text("report.css")
    linker_filtering = _load_static_image("linker_filtering.png")
    pareto_good = _load_static_image("pareto_good.png")
    pareto_one = _load_static_image("pareto_bad_one.png")
    pareto_many = _load_static_image("pareto_bad_many.png")
    barcode_qc_uri = _load_static_image("barcode_qc.png")
    low_lanes_correctable = _load_static_image("low_lanes_correctable.png")
    low_lanes_biological = _load_static_image("low_lanes_biological.png")
    high_lanes_correctable = _load_static_image("high_lanes_correctable.png")
    lane_failure = _load_static_image("lane_failure.png")

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
  {% if css_text %}
  <style>
{{ css_text }}
  </style>
  {% endif %}
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

          <details class="spoiler">
            <summary>Explain Linker Filtering</summary>
            <div class="spoiler-body">
              <p>
                Processing FASTQ files generated from a DBiT-seq experiment begins with two rounds of read-level filtering designed to remove sequences with low base quality or incorrect molecular structure.
                <br><br>
                This filtering step leverages the sequence identity of Ligation Linker 1 (L1) and Ligation Linker 2 (L2) (see figure below). These linker sequences are conserved regions that should be identical across all properly constructed molecules.
                <br><br>
                {% if linker_filtering %}
                  <div class="linker-example">
                    <img src="{{ linker_filtering }}" alt="linker_filtering" />
                  </div>
                {% endif %}
                <br>
                Because L1 and L2 are invariant and located at defined positions within each read, their sequence fidelity serves as an internal structural checkpoint. Reads with mismatches, indels, or truncated linker sequences are likely the result of sequencing errors, incomplete ligation, or improper library construction and are therefore removed.
                <br><br>
                Reads with 3 mismatches or more in either ligation linker sequence are removed from processing. Either filtering step removing >30% of the reads results in a CAUTION status. We recommend that the experiment proceed through processing--but it may be a canidate for re-sequencing if the downstream metrics do not meet desired criteria.
                </p>
                <br><br
                <p>From this filtering step, we compute the following metrics:</p>
                <ul>
                  <li>
                    <strong>Total Reads:</strong> The total number of input reads after subsampling and prior to filtering. This value should match the subsampling target unless the number of available raw reads was lower than the target.
                  </li>
                  <br>
                  <li>
                    <strong>Total with Linker:</strong> The number of reads containing the ligation linker sequence with fewer than three mismatches.
                  </li>
                  <br>
                  <li>
                    <strong>Percent Pass Filter:</strong> The percentage of input reads retained after filtering (Total with Linker ÷ Total Reads).
                  </li>
                  </li>
                  <br>
                  <li>
                    <strong>Number of Unique Barcodes:</strong>The number of distinct 8-mer sequences detected downstream of the linker. Because random sequencing errors introduce spurious variants, this value is typically higher than the expected number of barcodes. In well-performing experiments, however, unexpected 8-mers account for only a small fraction of reads and can generally be ignored (see Barcode Check section). discussion).
                  </li>
                  <br>
                  <li>
                    <strong>Number Barcodes with 90% of reads:</strong> Reads are grouped by barcode, and this metric reports how many barcodes collectively account for 90% of total reads. This value determines how many 8-mers are displayed in the Pareto plot (see Barcode Check below) and should typically be close to the expected number of barcodes.
                  </li>
                  <br>
                  <li>
                    <strong>Percent of Reads in Expected Barcodes:</strong> The percentage of total reads assigned to the barcode whitelist specified by the <code>barcode_set</code> parameter. This value is typically greater than 80%. See the Barcode Check section for additional details.
                  </li>
                </ul>
            </div>
          </details>
          <br>

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
          <details class="spoiler">
            <summary>Explain Barcode Check</summary>
              <div class="spoiler-body">
                <p>
                  This section screens for failure modes resulting from a mismatch between the molecular barcodes used in the DBiT-seq experiment and the barcodes selected for computational processing. It also checks for the presence of contaminating 8-mer sequences at the barcode location that may result from poor sequencing quality.
                </p>
                <br>
                <p>
                  We define the Barcode A and Barcode B 8-mer sequences as the 8 base pairs upstream of Ligation Linker 1 and Ligation Linker 2, respectively (only sequences with fewer than 3 mismatches are considered; see Filtering above). For each barcode (A and B), we create a table of read counts for each 8-mer sequence. We calculate the relative frequency of each barcode, sort by frequency, and generate the Pareto plot below showing the sequences that account for 90% of all reads.
                </p>
                <br>
                <p>
                  In the figure below, the left y-axis displays the frequency of each barcode, and the right y-axis shows the cumulative sum. Each point on the graph represents a detected 8-mer sequence. Expected 8-mers are colored green; the x-axis label for each point is comprised of the 8-mer sequence proceeded by either the channel for expected barcodes or 'nan' for unexpected.
                </p>
              <br>
              <p>
                Below, we present examples of common scenarios to aid interpretation.
              </p>
              <br>

              <h4>Success</h4>
                <br>
                <p>
                  This figure illustrates a successful experiment. The frequency plot (left axis) shows a smooth, steadily decreasing trend. The top-ranked 220 8-mers are all expected sequences (shown in green). Beyond these dominant barcodes, the frequency drops sharply, and the remaining 8-mers are unexpected (shown in black) and represent only a small fraction of total reads. The plot includes 8-mers cumulatively until 90% of input reads are accounted for.
                </p>
                <br>
                {% if pareto_good %}
                  <div class="lane-qc-example">
                    <img src="{{ pareto_good }}" alt="pareto_good" />
                  </div>
                {% endif %}

              <br><br>
              <h4>Wrong barcode_set (multiple barcodes)</h4>
                <br>
                <p>
                  The example below shows a run processed with an incorrect barcode set. Note the interspersed green and black points in the frequency plot, indicating a mixture of expected and unexpected barcode sequences. This can occur if, for example, the experiment used the bc220 barcode set but bc96 was selected during processing. If you observe this pattern in your run, verify that the correct barcode files were used.
                </p>
                <br>
                {% if pareto_many %}
                  <div class="lane-qc-example">
                    <img src="{{ pareto_many }}" alt="pareto_many" />
                  </div>
                {% endif %}
              <br><br>

              <h4>Wrong barcode_set (single barcode)</h4>
                <br>
                <p>
                  The example below shows a run with a single barcode mismatch. This typically occurs when the selected barcode set differs from the correct set by only one barcode—for example, using bc220 instead of bc220_05-OCT. If you observe this pattern in your run, verify that the correct barcode files were used.
                </p>
                <br>
                {% if pareto_one %}
                  <div class="lane-qc-example">
                    <img src="{{ pareto_one }}" alt="pareto_one" />
                  </div>
                {% endif %}
              <br>

            </div>
          </details>
          <br>

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
          <h2>Lane QC</h2>
          <details class="spoiler">
            <summary>Explain Lane QC</summary>
            <div class="spoiler-body">
              <p>
                In this section, we display the relative read counts for each
                barcode as a bar plot, organized by row and column. Expected
                barcodes—specified by the barcode_set parameter—are arranged
                left to right along the x-axis and separated into two groups:
                rows and columns (typically Barcode A and Barcode B,
                respectively). Column (X-coordinate) barcodes are ordered left
                to right, while row (Y-coordinate) barcodes are ordered top to
                bottom. Data for this figure can be found in the 'tables'
                directory of the outputs folder.  L1_counts.csv contains data
                for Barcode A, L2_counts.csv for Barcode B.
              </p>
              <p>
                The y-axis represents the fraction of total reads. For each
                lane, read counts across all pixels are summed and divided by
                the total read count, yielding the Fraction of Reads
                (see figure below).
              </p>
              {% if barcode_qc_uri %}
                <div class="lane-qc-figure">
                  <img src="{{ barcode_qc_uri }}" alt="barcode_qc" />
                </div>
              {% endif %}
              <br>
              <p>
                This QC plot serves as a quick diagnostic for identifying
                anomalous lanes. Barcodes with a read fraction greater than 2×
                the mean are flagged as potential high artifacts, while those
                with a fraction less than 0.5× the mean are flagged as
                potential low artifacts. In the bar plot, these potential
                artifacts are highlighted in red.
              </p>
              <p>
                Barcodes with unusually high or low read fractions may result
                from microfluidic artifacts. In many cases, these can be
                mitigated through computational smoothing (see documentations
                for 'barcodeqc preview' and 'barcodeqc smooth'). However, such
                patterns may also reflect genuine biological or tissue
                features, such as missing tissue or tissue boundaries.
              </p>
              <p>
                Below, we present examples of common scenarios to aid
                interpretation.
              </p>
              <br>
              <h4>Low Lanes</h4>
              <br>

                <p>Example 1</p>
                <p>
                  In the figure below, we present a barcode QC plot containing multiple correctable low lanes. This run would trigger a CAUTION status in the summary table at the top of the page.
                  <br><br>
                  The row count bar plot (bottom) contains two lanes with fractions less than half of the row mean and substantially lower than their neighboring lanes (red arrows). These anomalous lanes likely represent chip defects resulting from microfluidic flow irregularities. Both can be corrected using computational smoothing.
                  <br><br>
                  The top panel (columns) shows low lanes at the left and right edges of the assayed area. These “edge effects” are common artifacts and can also be corrected with smoothing.
                </p>

                {% if low_lanes_correctable %}
                  <div class="lane-qc-example">
                    <img src="{{ low_lanes_correctable }}" alt="low_lanes_correctable" />
                  </div>
                {% endif %}

                <br>

                <p>Example 2</p>
                <p>
                  In this example, the run contains “low lanes” that reflect genuine tissue features. These lanes would trigger a CAUTION status in the summary table at the top of the page, but they do not require corrective action.
                  <br><br>
                  The lanes in question are circled in green in the bottom bar plot and are nearly all below the 0.5× mean threshold. The right panel shows a spatial heatmap of the assay region colored by the log of read counts per pixel. The region corresponding to the low lanes is highlighted in green and contains either no tissue (left-most) or minimal tissue at the boundary.
                  <br><br>
                  Because low read counts are expected in these regions, the CAUTION status can be safely ignored.
                </p>

                {% if low_lanes_biological %}
                  <div class="lane-qc-example">
                    <img src="{{ low_lanes_biological }}" alt="low_lanes_biological" />
                  </div>
                {% endif %}

                <br>

              <h4>High Lanes</h4>
              <br>

                <p>Example 1</p>
                <p>
                  Here we show a run with four high lanes that can be corrected using computational smoothing. In the barcode plot, these lanes have a higher fraction of reads compared to their neighbors. This pattern likely represents a spatial artifact that can be resolved through smoothing and is expected to have minimal impact on downstream analysis.
                </p>
                <br>
                {% if high_lanes_correctable %}
                  <div class="lane-qc-example">
                    <img src="{{ high_lanes_correctable }}" alt="high_lanes_correctable" />
                  </div>
                {% endif %}

                <br>

                <h4>Contact Support</h4>
                <br>

                <p>Example 1</p>
                <p>
                  Below is the lane QC bar plot for a run with severe high and low artifacts. These artifacts display an alternating high/low pattern that makes computational smoothing unreliable. It is unlikely that this run can be adequately corrected.
                  <br><br>
                  If your run resembles this example, please contact your AtlasXomics support scientist for assistance.
                </p>
                <br>
                {% if lane_failure %}
                  <div class="lane-qc-example">
                    <img src="{{ lane_failure }}" alt="lane_failure" />
                  </div>
                {% endif %}
                <br>

            </div>
          </details>
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
        css_text=css_text,
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
        linker_filtering=linker_filtering,
        pareto_good=pareto_good,
        pareto_many=pareto_many,
        pareto_one=pareto_one,
        barcode_qc_uri=barcode_qc_uri,
        low_lanes_correctable=low_lanes_correctable,
        low_lanes_biological=low_lanes_biological,
        high_lanes_correctable=high_lanes_correctable,
        lane_failure=lane_failure,
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
