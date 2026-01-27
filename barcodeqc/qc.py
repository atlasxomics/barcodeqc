import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import subprocess

from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
from typing import Literal, Optional

import barcodeqc.utils as utils

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s",
    level=logging.INFO
)


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
        logging.info(f"Output dir {output_dir} not detected, making...")
        output_dir.mkdir(parents=True, exist_ok=True)

    if not r2_path.exists():
        raise FileNotFoundError(f"fastq file path does not exist: {r2_path}")

    if tissue_position_file is not None:
        if not tissue_position_file.exists():
            raise FileNotFoundError(
                f"Could not find tissue_postion file: {tissue_position_file}"
            )
        else:
            logging.info(f"Using tissue_postions_file: {tissue_position_file}")
    else:
        tissue_position_file = utils.BARCODE_PATHS[barcode_set]["positions"]
        logging.info(
            f"No tissue positions supplied, using: {tissue_position_file.name}"
        )

    # Run subsampling with seqtk
    ds_path = output_dir / f"ds_{sample_reads}.fastq.gz"
    ds_cmd = f"seqtk sample -s {random_seed} {r2_path} {sample_reads} | gzip > {ds_path}"

    logging.info(f"Running subsample command:\n{ds_cmd}")
    subprocess.run(ds_cmd, shell=True, check=True)
    logging.info("Completed subsampling")

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

    logging.info(f"Starting dmuxL1:\n{dmuxL1}")
    subprocess.run(dmuxL1, shell=True, check=True)
    logging.info("Completed dmuxL1")

    logging.info(f"Starting dmuxL2:\n{dmuxL2}")
    subprocess.run(dmuxL2, shell=True, check=True)
    logging.info("Completed dmuxL2")

    # SpatialTable
    spatial_table = utils.make_spatial_table(
        wc_linker1, wc_linker2, tissue_position_file
    )
    spatial_table_path = output_dir / "spatialTable.csv"
    spatial_table.to_csv(spatial_table_path, index=False)

    # For each dmux, process wildcard outputs.
    wcList = [wc_linker1, wc_linker2]
    logList = [log_linker1, log_linker2]
    expList = ["L1", "L2"]
    countTableList = []
    maxToNinety = -1
    outLines = []
    for wc, logF, eL in zip(wcList, logList, expList):
        logging.info(f"wcFile: {wc}\\tlog: {logF}\\tsample{eL}")

        try:
            df = utils.load_wc_file(wc)
        except ValueError as e:
            logging.error(f"{e}")
            exit(0)

        uniqueCounts = df['8mer'].value_counts()
        countTable = pd.DataFrame(uniqueCounts)
        countTable['frac_count'] = uniqueCounts / uniqueCounts.sum()

        countTable = countTable.sort_values(by=['frac_count'], ascending=False)
        countTable['cumulative_sum'] = countTable['frac_count'].cumsum()
        numToNinety = (countTable['cumulative_sum'] <= 0.9).sum()
        try:
            pctFor50 = f"{countTable.iloc[:50, 1].sum():.1%}"
            pctFor96 = f"{countTable.iloc[:96, 1].sum():.1%}"
        except:
            pctFor50 = 'NaN'
            pctFor96 = 'NaN'

        # mark bcs that are expected from whitelist
        intSet = (set(countTable.index.tolist()) & set(whitelist))

        # add column and set matches in set to True
        countTable['expectMer'] = False
        countTable['channel'] = -1
        countTable['merLabel'] = countTable.index

        for mer in list(intSet):
            countTable.loc[mer, 'expectMer'] = True

        countTableList.append(countTable)
        countTable.to_csv(output_dir / f'{eL}_counts.csv', index=False)
        totalReadFromExpected = countTable['frac_count'][countTable['expectMer']].sum()
        total_reads, adapter_reads = utils.parse_read_log(logF)

        out = f'\n######### Info for {os.path.basename(wc)}\n'
        out = out + f'Total Reads: {total_reads}  Reads with Adapter: {adapter_reads}'
        out = out + f'\nThe number of unique strings in the 8mer column is {len(uniqueCounts)}.\nNinety percent (90%) of the reads come from total of {numToNinety} 8mers.'
        out = out + f"\nTotal of {len(intSet)} out of {len(whitelist)} expected 8-mers accounted for {totalReadFromExpected:.1%} of the reads"
        out = out + f"\nTop 50 8mers represent {pctFor50} fraction of reads\nTop 96 8mers represent {pctFor96} fraction of reads"

        outLines.append(out)
        logging.info(out)

        if numToNinety > maxToNinety:
            maxToNinety = numToNinety

    # # For each of cutadapt result (L1, L2), identify hi/lo mers, save
    # # table as png if present, make barcode qc plot with hi/lo thresholds, ####
    # # and make barcode plate heatmap as png.
    # picPaths = []
    # tablePicPaths = []
    # for ctb, eL in zip(countTableList, expList):
    #     expectedTable = ctb.copy()
    #     expectedTable['8mer'] = expectedTable.index
    #     expectedTable = expectedTable[
    #         expectedTable.channel != -1
    #     ].sort_values(by=['channel'], ascending=True)
    #     expectedTable['hiWarn'] = False
    #     expectedTable['loWarn'] = False

    #     meanCount = expectedTable['frac_count'].mean()
    #     upperCut = 2 * meanCount
    #     lowerCut = 0.5 * meanCount

    #     expectedTable.loc[expectedTable['frac_count'] > upperCut, 'hiWarn'] = True
    #     expectedTable.loc[expectedTable['frac_count'] < lowerCut, 'loWarn'] = True

    #     totalHiWarn = expectedTable['hiWarn'].sum()
    #     totalLoWarn = expectedTable['loWarn'].sum()
    #     totalMers = len(expectedTable['frac_count'])
    #     out = f'\n######### Info for {eL}\n'
    #     out = out+f'Total hiWarn: {totalHiWarn}\tTotal loWarn: {totalLoWarn}'
    #     out = out+f'\nPct hiWarn: {(totalHiWarn/totalMers):.3f}\tPct loWarn: {(totalLoWarn/totalMers):.3f}'

    #     outLines.append(out)
    #     logging.info(out)

    #     # save hi/loWarn table
    #     if (totalHiWarn + totalLoWarn) > 0:  # Only export if there are hi or lows

    #         # Create a subset of the expectedTable where hiWarn and loWarn are true
    #         subset_expectedTable = expectedTable.loc[expectedTable['hiWarn'] | expectedTable['loWarn']]
    #         subset_expectedTable.to_csv(
    #             f'{output_dir}{eL}_hiLoWarn.csv', index=False
    #         )

    #         # After generating the hi/loWarn table, save png of table
    #         hi_lo_warn_table_path = f'{output_dir}{eL}_hiLoWarn.png'
    #         subset_expectedTable.loc[subset_expectedTable['hiWarn'] == True, 'hiOrLo'] = 'high'
    #         subset_expectedTable.loc[subset_expectedTable['loWarn'] == True, 'hiOrLo'] = 'low'
    #         expectColList = ['merLabel', 'plate', 'plateCol', 'plateRow', 'count', 'hiOrLo', 'fracCount']
    #         utils.save_table_as_image(
    #             subset_expectedTable[expectColList],
    #             hi_lo_warn_table_path
    #         )

    #         logging.info(f"Saved hi/loWarn table as image: {hi_lo_warn_table_path}")
    #         tablePicPaths.append(hi_lo_warn_table_path)

    #     # make bar plot (barcode QC plot)
    #     fig, ax = plt.subplots(figsize=(20, 5))
    #     expectedTable.plot(
    #         kind='bar',
    #         x='channel',
    #         y='frac_count',
    #         ax=ax,
    #         logy=False,
    #         legend=False
    #     )
    #     ax.set_xticks(expectedTable.channel)
    #     ax.set_xticklabels(expectedTable.merLabel, rotation=90)
    #     ax.axhline(y=meanCount, color='tab:gray', linestyle='--', linewidth=0.5, label=f'mean: {meanCount:.3f}')
    #     ax.axhline(y=meanCount/2, color='tab:red', linestyle='--', linewidth=0.5, label=f'half mean: {lowerCut:.3f}')
    #     ax.axhline(y=meanCount*2, color='tab:red', linestyle='--', linewidth=0.5, label=f'dbl mean: {upperCut:.3f}')
    #     fig.legend(loc="upper right")
    #     fig.tight_layout()

    #     # saving plot
    #     barplot_name = f'{output_dir}{eL}_plot.png'
    #     logging.info(f"Saving plot {barplot_name}...")
    #     fig.savefig(barplot_name)
    #     picPaths.append(barplot_name)
    #     logging.info(f"Completed {barplot_name}")

    # # Do some on/off tissue calcs ####
    # totalTix = len(spatial_table['on_off'])
    # totalOnTis = len(spatial_table.loc[spatial_table['on_off'] == 1]['on_off'])
    # totalOffTis = len(spatial_table.loc[spatial_table['on_off'] == 0]['on_off'])
    # countsOnTis = spatial_table.loc[spatial_table['on_off'] == 1]['count'].sum()
    # countsOffTis = spatial_table.loc[spatial_table['on_off'] == 0]['count'].sum()
    # countsPerTixOnTis = countsOnTis / totalOnTis

    # # generate density plot for on/off tissue pixels
    # density_plotPath = 'output_dir / 'denseon_off.png'

    # utils.create_density_plot(
    #     spatial_table,
    #     density_plotPath,
    #     'count',
    #     'on_off',
    #     log10=True,
    #     x_label='total counts',
    #     y_label='density'
    # )
    # picPaths.append(density_plotPath)

    # try:
    #     countsPerTixOffTis = countsOffTis / totalOffTis
    # except ZeroDivisionError:
    #     countsPerTixOffTis = 0

    # fractCPToffTiss = countsPerTixOffTis / countsPerTixOnTis
    # ratioOffvOn = countsOffTis / countsOnTis

    # logging.info(
    #     f"TixelStats: {totalTix=}, {totalOnTis=}, {totalOffTis=}, \
    #     {countsOnTis=}, {countsOffTis=}, {ratioOffvOn=}, \
    #     {countsPerTixOnTis=}, {countsPerTixOffTis=}, {fractCPToffTiss=}"
    # )

    # valList = [
    #     totalTix, totalOnTis, totalOffTis, countsOnTis,
    #     countsOffTis, countsPerTixOnTis, countsPerTixOffTis,
    #     fractCPToffTiss, ratioOffvOn
    # ]

    # nameList = ['totalTix', 'totalOnTis', 'totalOffTis', 'countsOnTis',
    #             'countsOffTis', 'countsPerTixOnTis', 'countsPerTixOffTis',
    #             'fractCPToffTiss', 'ratioOffvOn']

    # mTable = utils.variables_to_dataframe(valList, nameList)
    # mTablePath = output_dir / "on_offTissueTable_mqc.csv"
    # mTable.to_csv(mTablePath, index=False)

    # # write html with plate/plot figures, append table pics to the end of ####
    # #  picPaths ####
    # picPaths = picPaths + tablePicPaths
    # reportNoteHTML = f"<br><b>Reads with L1:</b><br>{tRead}<br><b>Reads with L1&L2:</b><br>{aRead}"
    # utils.generate_html_with_embedded_images(
    #     picPaths, output_dir, sample_name, noteHTML=reportNoteHTML
    # )

    # # Make pareto chart of barcode abundances, save as _output.pdf ####
    # with PdfPages(output_dir / 'output.pdf') as pdf:
    #     for tb, wc, ol in zip(countTableList, wcList, outLines):
    #         # define subplots
    #         startI = 0
    #         endI = maxToNinety + int(maxToNinety / 2)
    #         fig, ax = plt.subplots(figsize=(20, 10))

    #         ax.plot(tb['fracCount'][startI:endI], marker='x',
    #                 linewidth=1, label='Fraction total 8mer')

    #         # define second y-axis that shares x-axis with current plot
    #         ax2 = ax.twinx()

    #         # add second line to plot
    #         ax2.plot(
    #             tb['cumulative_sum'][startI:endI], marker='.', linewidth=1,
    #             color='tab:orange', label='Cumlative Fraction'
    #         )

    #         # add second y-axis label
    #         labels0 = list(tb.index[startI:endI])
    #         labels2 = list(tb.index[startI:endI] + "_" + tb.channel[startI:endI].astype(str).str.zfill(2))

    #         chanLabels = list("_" + (tb.channel[startI:endI] + 1).astype(str).str.zfill(2))
    #         labels3 = [f"{lb}{cb}" for lb, cb in zip(labels0, chanLabels)]

    #         labels = labels3
    #         ax.set_xticks(labels0)
    #         ax.set_xticklabels(labels, rotation=90)

    #         s = tb.expectMer[startI:endI]
    #         makeGreen = np.where(s)[0].tolist()
    #         for mg in makeGreen:
    #             ax.get_xticklabels()[mg].set_color("green")

    #         ax2.axhline(y=.9, color='tab:gray', linestyle='--', linewidth=0.5, label='90 percent')
    #         fig.legend(loc="center right")
    #         plt.title(f"plot: {wc}")
    #         plt.annotate(ol, xy=(0.6, 0.6), xycoords='axes fraction', ha='left', va='center', ma='left')
    #         fig.subplots_adjust(bottom=0.2)
    #         plt.close()

    #         pdf.savefig(fig)

    return Path(spatial_table_path)
