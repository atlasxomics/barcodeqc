from __future__ import annotations

import logging
import pandas as pd
import subprocess

from pathlib import Path
from typing import Literal, Optional

import barcodeqc.files as files
import barcodeqc.plots as plots
import barcodeqc.utils as utils

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

    # Check fastq exists, if formatted propery
    # ensure sample_reads is a reasonable number

    output_dir = Path.cwd() / sample_name
    if not output_dir.exists():
        logger.debug(f"Output dir {output_dir} not detected, making...")
        output_dir.mkdir(parents=True, exist_ok=True)

    if not r2_path.exists():
        raise FileNotFoundError(f"fastq file path does not exist: {r2_path}")

    if tissue_position_file is not None:
        if not tissue_position_file.exists():
            raise FileNotFoundError(
                f"Could not find tissue_postion file: {tissue_position_file}"
            )
        else:
            logger.debug(f"Using tissue_postions_file: {tissue_position_file}")
    else:
        tissue_position_file = utils.BARCODE_PATHS[barcode_set]["positions"]
        logger.debug(
            f"No tissue positions supplied, using: {tissue_position_file.name}"
        )

    bca_file = utils.BARCODE_PATHS[barcode_set]["bca"]
    bca_positions = files.open_barcode_file(bca_file)

    bcb_file = utils.BARCODE_PATHS[barcode_set]["bcb"]
    bcb_positions = files.open_barcode_file(bcb_file)

    # Run subsampling with seqtk
    ds_path = output_dir / f"ds_{sample_reads}.fastq.gz"
    ds_cmd = f"seqtk sample -s {random_seed} {r2_path} {sample_reads} | gzip > {ds_path}"

    logger.info("Running subsample ")
    logger.debug(ds_cmd)
    subprocess.run(ds_cmd, shell=True, check=True)
    logger.info("Completed subsampling")

    # Run filtering with cutadapt
    linker1 = "NNNNNNNNGTGGCCGATGTTTCGCATCGGCGTACGACT"
    linker2 = "NNNNNNNNATCCACGTGCTTGAGAGGCCAGAGCATTCG"

    wc_linker1 = output_dir / "cutadapt_wc_L1.txt"
    log_linker1 = output_dir / "cutadapt_L1.log"
    dmuxL1 = (
        f"cutadapt "
        f"-g linker1={linker1} "
        f"-o /dev/null "
        "--action=lowercase --cores 30 "
        "--no-indels -e 5 "
        f"--wildcard-file {wc_linker1} "
        f"{ds_path} "
        f"> {log_linker1}"
    )

    wc_linker2 = output_dir / "cutadapt_wc_L2.txt"
    log_linker2 = output_dir / "cutadapt_L2.log"
    dmuxL2 = (
        f"cutadapt "
        f"-g linker2={linker2} "
        f"-o /dev/null "
        "--action=lowercase --cores 30 "
        "--no-indels -e 5 "
        f"--wildcard-file {wc_linker2} "
        f"{ds_path} "
        f"> {log_linker2}"
    )

    logger.info("Starting dmuxL1")
    logger.debug(dmuxL1)
    subprocess.run(dmuxL1, shell=True, check=True)
    logger.info("Completed dmuxL1")

    logger.info("Starting dmuxL2")
    logger.debug(dmuxL2)
    subprocess.run(dmuxL2, shell=True, check=True)
    logger.info("Completed dmuxL2")

    # SpatialTable
    spatial_table = utils.make_spatial_table(
        wc_linker1, wc_linker2, tissue_position_file
    )
    spatial_table_path = output_dir / "spatialTable.csv"
    spatial_table.to_csv(spatial_table_path, index=False)

    # For each dmux, process wildcard outputs.
    wc_list = [wc_linker1, wc_linker2]
    logList = [log_linker1, log_linker2]
    expList = ["L1", "L2"]
    bc_list = [bca_positions, bcb_positions]
    row_col = ["row", "col"]
    countTableList = []
    pic_paths = []
    maxToNinety = -1

    for wc, logF, eL, bcl, rc in zip(wc_list, logList, expList, bc_list, row_col):

        logger.info(f"Processing 8mer counts for {eL}")
        logger.debug(f"wcFile: {wc}\\tlog: {logF}\\tsample: {eL}")

        whitelist = bcl["sequence"]
        wc_df = files.load_wc_file(wc)

        unique_counts = wc_df['8mer'].value_counts()
        count_table = pd.DataFrame(unique_counts)
        count_table['frac_count'] = unique_counts / unique_counts.sum()

        count_table = count_table.sort_values(by=['frac_count'], ascending=False)
        count_table['cumulative_sum'] = count_table['frac_count'].cumsum()
        numToNinety = (count_table['cumulative_sum'] <= 0.9).sum()

        try:
            pctFor50 = f"{count_table.iloc[:50, 1].sum():.1%}"
            pctFor96 = f"{count_table.iloc[:96, 1].sum():.1%}"
        except:
            pctFor50 = 'NaN'
            pctFor96 = 'NaN'

        # mark bcs that are expected from whitelist
        expected_bcs = (set(count_table.index.tolist()) & set(whitelist))

        # Merge row/col data from bcl table
        count_table = count_table.join(
            bcl.set_index('sequence')[['row', 'col']],
            how='left'
        )
        count_table['sequence'] = count_table.index

        count_table['channel'] = -1
        # populate channel with specified row_col value
        count_table['channel'] = count_table[rc]

        count_table['expectMer'] = count_table['sequence'].isin(expected_bcs)

        count_table.to_csv(output_dir / f'{eL}_counts.csv', index=True)
        total_read_from_expected = count_table['frac_count'][count_table['expectMer']].sum()

        countTableList.append(count_table)

        total_reads, adapter_reads = utils.parse_read_log(logF)
        out = f'\n######### Info for {wc.name}\n'
        out = out + f'Total Reads: {total_reads}  Reads with Adapter: {adapter_reads}'
        out = out + f'\nThe number of unique strings in the 8mer column is {len(unique_counts)}.\nNinety percent (90%) of the reads come from total of {numToNinety} 8mers.'
        out = out + f"\nTotal of {len(expected_bcs)} out of {len(whitelist)} expected 8-mers accounted for {total_read_from_expected:.1%} of the reads"
        out = out + f"\nTop 50 8mers represent {pctFor50} fraction of reads\nTop 96 8mers represent {pctFor96} fraction of reads"

        logger.debug(out)

        if numToNinety > maxToNinety:
            maxToNinety = numToNinety

        logger.info(f"Identifying hi/lo barcodes for {eL}")
        bc_table = count_table.copy()

        bc_table = bc_table[bc_table["expectMer"]].sort_values(
            by=['channel'], ascending=True
        )

        bc_table['hiWarn'] = False
        bc_table['loWarn'] = False

        mean = bc_table['frac_count'].mean()
        upperCut = 2 * mean
        lowerCut = 0.5 * mean

        bc_table.loc[bc_table['frac_count'] > upperCut, 'hiWarn'] = True
        bc_table.loc[bc_table['frac_count'] < lowerCut, 'loWarn'] = True

        totalHiWarn = bc_table['hiWarn'].sum()
        totalLoWarn = bc_table['loWarn'].sum()
        totalMers = len(bc_table)

        out = f'\n######### Info for {eL}\n'
        out = out+f'Total hiWarn: {totalHiWarn}\tTotal loWarn: {totalLoWarn}'
        if totalMers > 0:
            out = out + f'\nPct hiWarn: {(totalHiWarn / totalMers):.3f}\tPct loWarn: {(totalLoWarn / totalMers):.3f}'
        else:
            out = out + '\nPct hiWarn: N/A\tPct loWarn: N/A'

        logger.debug(out)

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
        logger.info("Barplot saved...")

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
        logger.info("Pareto plot saved...")

    # Do some on/off tissue calcs
    logger.info("Calculating on/off tissue stats...")
    on_df = spatial_table.loc[spatial_table['on_off'] == 1]
    off_df = spatial_table.loc[spatial_table['on_off'] == 0]

    total_pix = len(spatial_table)
    total_on = len(on_df)
    total_off = len(off_df)
    counts_on = on_df['count'].sum()
    counts_off = off_df['count'].sum()
    counts_per_pix_on = counts_on / total_on
    if total_off > 0:
        counts_per_pix_off = counts_off / total_off
    else:
        counts_per_pix_off = 0
    frac_per_pix_off_on = counts_per_pix_off / counts_per_pix_on
    ratio_off_on = counts_off / counts_on

    logger.debug(
        f"Pixel Stats: {total_pix=}, {total_on=}, {total_off=}, "
        f"{counts_on=}, {counts_off=}, {ratio_off_on=}, "
        f"{counts_per_pix_on=}, {counts_per_pix_off=}, {frac_per_pix_off_on=}"
    )

    onoff_dict = {
        "metric": ["total_pix", "total_on", "total_off", "counts_on", "counts_off", "ratio_off_on", "counts_per_pix_on", "counts_per_pix_off", "frac_per_pix_off_on"],
        "value":  [total_pix, total_on, total_off, counts_on, counts_off, ratio_off_on, counts_per_pix_on, counts_per_pix_off, frac_per_pix_off_on]
    }
    onoff_df = pd.DataFrame(onoff_dict)

    onoff_path = output_dir / "onoff_tissue_table.csv"
    onoff_df.to_csv(onoff_path, index=False)

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
    # write html with plate/plot figures
    report_note = f"<br><b>Reads with L1:</b><br>{total_reads}<br><b>Reads with L1&L2:</b><br>{adapter_reads}"
    utils.generate_html_with_embedded_images(
        pic_paths, output_dir, sample_name, noteHTML=report_note
    )
    logger.info("html report finished.")

    return Path(spatial_table_path)
