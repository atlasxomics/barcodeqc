from __future__ import annotations

import logging
from pathlib import Path
from typing import Literal, Optional

import barcodeqc.files as files
import barcodeqc.plots as plots
import barcodeqc.utils as utils
from barcodeqc.logging_helpers import (
    format_hilo_metrics, format_wildcard_metrics
)
from barcodeqc.qc_config import QCConfig
from barcodeqc.qc_steps import (
    build_count_table,
    build_spatial_table,
    compute_hi_lo_qc,
    compute_onoff_metrics,
    write_onoff_table,
)
from barcodeqc.qc_steps import ensure_output_dir, run_cutadapt, run_subsample

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

    logger.debug(
        f"Using tissue_postions_file: {config.tissue_position_file}"
    )

    bca_file = utils.BARCODE_PATHS[barcode_set]["bca"]
    bca_positions = files.open_barcode_file(bca_file)

    bcb_file = utils.BARCODE_PATHS[barcode_set]["bcb"]
    bcb_positions = files.open_barcode_file(bcb_file)

    ds_path = run_subsample(config)

    wc_linker1, log_linker1, wc_linker2, log_linker2 = run_cutadapt(
        ds_path, config.output_dir
    )

    spatial_table, spatial_table_path = build_spatial_table(
        wc_linker1,
        wc_linker2,
        config.tissue_position_file,
        config.output_dir,
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
            pctFor50,
            pctFor96,
            whitelist
        ) = build_count_table(wc, bcl, rc)

        count_table.to_csv(output_dir / f'{eL}_counts.csv', index=True)
        total_read_from_expected = count_table['frac_count'][
            count_table['expectMer']
        ].sum()

        countTableList.append(count_table)

        total_reads, adapter_reads = utils.parse_read_log(logF)
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
                pctFor50,
                pctFor96,
            )
        )

        if numToNinety > maxToNinety:
            maxToNinety = numToNinety

        logger.info(f"Identifying hi/lo barcodes for {eL}...")
        bc_table, totalHiWarn, totalLoWarn, totalMers = compute_hi_lo_qc(
            count_table
        )

        logger.debug(
            format_hilo_metrics(eL, totalHiWarn, totalLoWarn, totalMers)
        )

        # Only export if there are hi/lows
        if (totalHiWarn + totalLoWarn) > 0:

            subset_expectedTable = bc_table.loc[
                bc_table['hiWarn'] | bc_table['loWarn']
            ]
            subset_expectedTable.to_csv(
                output_dir / f"{eL}_hiLoWarn.csv", index=False
            )

        logger.info(f"Saving barcode barplot for {eL}...")
        barplot_path = plots.hilo_plot(
            bc_table,
            "channel",
            "frac_count",
            "sequence",
            output_dir,
            f"{eL}_barplot.png",
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
            output_dir,
            f"{eL}_pareto.png",
        )
        pic_paths.append(pareto_path)
        logger.info("Pareto plot saved.")

    # Do some on/off tissue calcs
    logger.info("Calculating on/off tissue stats...")
    onoff_df = compute_onoff_metrics(spatial_table)
    _ = write_onoff_table(onoff_df, output_dir)

    # generate density plot for on/off tissue pixels
    density_path = output_dir / "dense_on_off.png"
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
    utils.generate_html_with_embedded_images(
        pic_paths, output_dir, sample_name
    )
    logger.info("html report finished.")

    return Path(spatial_table_path)
