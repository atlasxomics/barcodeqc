from __future__ import annotations

import logging
import pandas as pd
from pathlib import Path
from typing import Literal, Optional

import barcodeqc.files as files
import barcodeqc.paths as paths
import barcodeqc.plots as plots
import barcodeqc.report as report
import barcodeqc.utils as utils
from barcodeqc.logging import (
    format_hilo_metrics, format_wildcard_metrics
)
from barcodeqc.config import QCConfig
from barcodeqc.steps import (
    build_count_table,
    build_spatial_table,
    barcode_check_status,
    compute_hi_lo_qc,
    compute_onoff_metrics,
    ensure_output_dir,
    lane_status,
    linker_conservation_status,
    run_cutadapt,
    run_subsample,
    write_onoff_table,
)

logger = logging.getLogger(__name__)


def qc(
    sample_name: str,
    r2_path: Path,
    barcode_set: Literal[
        "bc50", "bc96", "fg96", "bc220", "bc220_05-OCT", "bc220_20-MAY"
    ],
    sample_reads: int,
    random_seed: int,
    tissue_position_file: Optional[Path]
) -> Path:

    # Setup
    config = QCConfig.from_args(
        sample_name=sample_name,
        r2_path=r2_path,
        barcode_set=barcode_set,
        sample_reads=sample_reads,
        random_seed=random_seed,
        tissue_position_file=tissue_position_file,
    )
    config.validate()
    ensure_output_dir(config.output_dir)

    output_dir = config.output_dir
    figures_dir = output_dir / "figures"
    tables_dir = output_dir / "tables"
    logs_dir = output_dir / "logs"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    logger.debug(
        f"Using tissue_postions_file: {config.tissue_position_file}"
    )

    raw_reads = None
    try:
        raw_reads = utils.count_fastq_reads(config.r2_path)
    except Exception as exc:
        logger.warning(
            "Failed to count raw reads in %s: %s", config.r2_path, exc
        )

    bca_file = paths.BARCODE_PATHS[barcode_set]["bca"]
    bca_positions = files.open_barcode_file(bca_file)

    bcb_file = paths.BARCODE_PATHS[barcode_set]["bcb"]
    bcb_positions = files.open_barcode_file(bcb_file)

    # Run subsample command with seqtk
    ds_path = run_subsample(config, output_dir)

    # Filtering linkers and parse barcodes with cutadapt
    wc_linker1, log_linker1, wc_linker2, log_linker2 = run_cutadapt(
        ds_path, logs_dir
    )

    # Build spatial table
    spatial_table, spatial_table_path = build_spatial_table(
        wc_linker1,
        wc_linker2,
        config.tissue_position_file,
        tables_dir,
    )

    # For each dmux, process wildcard outputs.
    wc_list = [wc_linker1, wc_linker2]
    logList = [log_linker1, log_linker2]
    expList = ["L1", "L2"]
    bc_list = [bca_positions, bcb_positions]
    row_col = ["row", "col"]
    countTableList = []
    pic_paths = []
    maxToNinety = -1
    linker_metrics: dict[str, dict[str, str | int | float]] = {}
    linker_status: dict[str, tuple[str, float]] = {}
    barcode_status: dict[str, tuple[str, int]] = {}
    hi_lane_statuses: list[str] = []
    lo_lane_statuses: list[str] = []

    for wc, logF, eL, bcl, rc in zip(
        wc_list, logList, expList, bc_list, row_col
    ):

        logger.info(f"Processing 8mer counts for {eL}")
        logger.debug(f"wcFile: {wc}\tlog: {logF}\tsample: {eL}")

        (
            count_table,
            unique_counts,
            expected_bcs,
            numToNinety,
            whitelist
        ) = build_count_table(wc, bcl, rc)

        count_table.to_csv(tables_dir / f'{eL}_counts.csv', index=True)
        total_read_from_expected = count_table['frac_count'][
            count_table['expectMer']
        ].sum()

        countTableList.append(count_table)

        total_reads_str, adapter_reads_str = utils.parse_read_log(logF)
        total_reads = int(total_reads_str)
        adapter_reads = int(adapter_reads_str)
        logger.debug(
            format_wildcard_metrics(
                wc.name,
                total_reads,
                adapter_reads,
                len(unique_counts),
                numToNinety,
                len(expected_bcs),
                len(whitelist),
                total_read_from_expected,
            )
        )
        linker_metrics[eL] = {
            "Total Reads": total_reads,
            "Total with Linker": adapter_reads,
            "Percent Pass Filtering": (
                f"{(adapter_reads / total_reads):.1%}"
                if total_reads > 0
                else "NA"
            ),
            "Number of Unique Barcodes": len(unique_counts),
            "Number Barcodes with 90% of reads": numToNinety,
            "Percent reads in expected barcodes": f"{total_read_from_expected:.1%}",
        }
        linker_status[eL] = linker_conservation_status(
            total_reads,
            adapter_reads,
        )
        barcode_status[eL] = barcode_check_status(count_table)

        if numToNinety > maxToNinety:
            maxToNinety = numToNinety

        logger.info(f"Identifying hi/lo barcodes for {eL}...")
        bc_table, totalHiWarn, totalLoWarn, totalMers = compute_hi_lo_qc(
            count_table
        )

        logger.debug(
            format_hilo_metrics(eL, totalHiWarn, totalLoWarn, totalMers)
        )
        hi_lane_statuses.append(lane_status(bc_table, "hiWarn"))
        lo_lane_statuses.append(
            lane_status(bc_table, "loWarn", edge_adjacent_ok=True)
        )

        # Only export if there are hi/lows
        if (totalHiWarn + totalLoWarn) > 0:

            subset_expectedTable = bc_table.loc[
                bc_table['hiWarn'] | bc_table['loWarn']
            ]
            subset_expectedTable.to_csv(
                tables_dir / f"{eL}_hiLoWarn.csv", index=False
            )

        logger.info(f"Saving barcode barplot for {eL}...")
        barplot_path = plots.hilo_plot(
            bc_table,
            "channel",
            "frac_count",
            "sequence",
            figures_dir,
            f"{eL}_barplot.html",
        )
        pic_paths.append(barplot_path)
        logger.info("Barplot saved.")

        # Make pareto chart of barcode abundances, save as _output.pdf ####
        logger.info(f"Saving pareto plot for {eL}")
        pareto_path = plots.pareto_plot(
            count_table,
            "frac_count",
            "cumulative_sum",
            "expectMer",
            "channel",
            wc,
            maxToNinety,
            figures_dir,
            f"{eL}_pareto.html",
        )
        pic_paths.append(pareto_path)
        logger.info("Pareto plot saved.")

    # Do some on/off tissue calcs
    logger.info("Calculating on/off tissue stats...")
    onoff_df = compute_onoff_metrics(spatial_table)
    _ = write_onoff_table(onoff_df, tables_dir)

    # generate density plot for on/off tissue pixels
    density_path = figures_dir / "dense_on_off.html"
    plots.create_density_plot(
        spatial_table,
        density_path,
        "count",
        "on_off",
        log10=True,
        x_label="total counts",
        y_label="density"
    )
    pic_paths.append(density_path)

    logger.info("on/off tissue stats finished.")

    logger.info("Generating html report...")
    if "CONTACT SUPPORT" in hi_lane_statuses:
        hi_lane_summary = "CONTACT SUPPORT"
    elif "ACTION REQUIRED" in hi_lane_statuses:
        hi_lane_summary = "ACTION REQUIRED"
    else:
        hi_lane_summary = "PASS"

    if "CONTACT SUPPORT" in lo_lane_statuses:
        lo_lane_summary = "CONTACT SUPPORT"
    elif "ACTION REQUIRED" in lo_lane_statuses:
        lo_lane_summary = "ACTION REQUIRED"
    else:
        lo_lane_summary = "PASS"

    ratio_row = onoff_df.loc[onoff_df["metric"] == "ratio_off_on", "value"]
    off_tissue_ratio = (
        f"{ratio_row.iloc[0]:.3f}" if not ratio_row.empty else "NA"
    )

    summary_table = pd.DataFrame(
        {
            "metric": [
                "Linker 1 Filter",
                "Linker 2 Filter",
                "Barcode A Check",
                "Barcode B Check",
                "HIGH Lanes",
                "LOW Lanes",
                "Off-tissue Ratio",
            ],
            "status": [
                linker_status.get("L1", ("NA", 0.0))[0],
                linker_status.get("L2", ("NA", 0.0))[0],
                barcode_status.get("L1", ("NA", 0))[0],
                barcode_status.get("L2", ("NA", 0))[0],
                hi_lane_summary,
                lo_lane_summary,
                off_tissue_ratio,
            ],
            "description": [
                "Reads are filtered on the sequence identity of the first ligation linker (L1); reads with more than 3 mismatches are removed from processing. PASS: >70% of reads kept. A low passing rate can indicate poor quality sequencing.",
                "Reads are filtered on the sequence identity of the second ligation linker (L2); reads with more than 3 mismatches are removed from processing. PASS: >70% of reads kept. A low passing rate can indicate poor quality sequencing.",
                "Barcode A sequences are extracted and compared against the user-defined whitelist. PASS: No unexpected sequences in the top 100 sequences(sorted by read count); CAUTION: >1 unexpected sequences.  Unexpected sequences can indicate a mismatch between the barcodes used and the whitelist selected for processing.",
                "Barcode B sequences are extracted and compared against the user-defined whitelist. PASS: No unexpected sequences in the top 100 sequences (sorted by read count); CAUTION: >1 unexpected sequences.  Unexpected sequences can indicate a mismatch between the barcodes used and the whitelist selected for processing.",
                "Reads with >2x the mean read count per row/col are flagged: PASS: no high lanes; ACTION REQUIRED: one or more non-adjacent high lane(s); CONTACT SUPPORT: adjacent high lanes. High lanes can be remediated with computational smoothing.",
                "Reads with <0.5x the mean read count per row/col are flagged: PASS: no low lanes; ACTION REQUIRED: one or more non-adjacent low lane(s) OR adjacent low lanes on edge of assay area (please check for off-tissue); CONTACT SUPPORT: adjacent low lanes internal to assay area (not on edge). Non-adjacent low lanes can be remediated with computational filling.",
                "Ratio of reads from off-tissue pixels to on-tissue pixels. On/off pixels are defined by the user-supplied tissue_positions_file; if no file was supplied, all pixels are considered 'on-tissue' and the ratio is 0.  A high ratio can indicate incorrect on/off tissue assignment in AtlasXBrowser or procedural artifacts."
            ]
        }
    )
    report.print_summary_table(summary_table)
    input_params = [
        {"label": "Sample Name", "value": sample_name},
        {"label": "Barcode File", "value": barcode_set},
        {
            "label": "Raw Reads",
            "value": f"{raw_reads:,}" if raw_reads is not None else "NA",
        },
    ]
    report.write_input_params(input_params, tables_dir)
    report.generate_report(
        figure_paths=pic_paths,
        output_dir=output_dir,
        sample_name=sample_name,
        summary_table=summary_table,
        linker_metrics=linker_metrics,
        onoff_table=onoff_df,
        input_params=input_params,
        file_tag="bcQC",
        table_dir=tables_dir,
    )
    logger.info("html report finished.")

    return Path(spatial_table_path)
