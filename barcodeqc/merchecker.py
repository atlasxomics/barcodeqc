import argparse
import json
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import subprocess

from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

from barcodeqc.bc_sets import bc_sets
import barcodeqc.utils as utils

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s",
    level=logging.INFO
)

repo_root = os.path.abspath(os.path.dirname(__file__))

output_dir = Path.cwd()
scripts_dir = os.path.join(repo_root, "barcodeqc")

parser = argparse.ArgumentParser(description="barcodeqc")

parser.add_argument(
    "-b",
    "--bcset",
    required=True,
    help="Barcode Set: bc96|bc50|bc210|fg9|fg210p30"
)
parser.add_argument(
    "-f",
    "--fq_path",
    required=True,
    help="path to Read2 FASTQ.GZ file"
)
parser.add_argument(
    "-r",
    "--sample_reads",
    required=False,
    default=1000000,
    help="downsampleReads default=10Million"
)
parser.add_argument(
    "-s",
    "--random_seed",
    required=False,
    default=100,
    help="Seed for randomization"
)
parser.add_argument(
    "-n",
    "--sample_name",
    required=True,
    default="{}",
)
parser.add_argument(
    "-t",
    "--tissue_position_file",
    required=False,
    default="none",
)

args = vars(parser.parse_args())
logging.info(f"args: {args}")

fastq_path = args['fq_path']
outReads = args['sample_reads']
seed = args['random_seed']
bcSet = args['bcset']
expName = args['sample_name']
tissue_pos_path = args['tissue_position_file']

# set barcode specific params and files
cfg = bc_sets.get(bcSet)
if cfg is None:
    logging.warning(f"No valid barcode set specified [bcSet: {bcSet}]")
    cfg = bc_sets['bc220']  # default to bc220
    logging.info("Defaulting to bc220 settings")
merPath = cfg["merPath_template"].format(merDir={scripts_dir})
plateList = cfg["plateList"]
logging.info(cfg["title"].format(bcSet=bcSet))

merTable = pd.read_csv(merPath, sep='\t')

isExist = os.path.exists(fastq_path)
if (not isExist):
    logging.error(f"FastQ file path does not exist: {fastq_path}")
else:
    logging.info(f"FastQ file path: {fastq_path}")

try:
    logging.info(f"ExpName: {expName}")
except:
    expName = os.path.basename(fastq_path).split('.')[0]
    logging.info(f"ExpName (auto): {expName}")

expDir = f"{output_dir}{expName}/"
subprocess.run(f'mkdir {expDir}', shell=True)

# Save the args as a JSON file
with open(f"{expDir}args.json", 'w') as json_file:
    json.dump(args, json_file, indent=4)

# Setup commands for seqtk, cutadapt ####
ds_fq = f"{fastq_path}_ds{outReads}.fq.gz"

linker1 = "NNNNNNNNGTGGCCGATGTTTCGCATCGGCGTACGACT"
linker2 = "NNNNNNNNATCCACGTGCTTGAGAGGCCAGAGCATTCG"

ds_cmd = f"seqtk sample -s{seed} {fastq_path} {outReads} | gzip > {ds_fq}"

dmuxL1 = (
    f"cutadapt "
    f"-g linker1={linker1} "
    f"-o \"{expDir}found-{expName}-{{name}}.fastq.gz\" "
    "--action=lowercase --cores 30 "
    "--no-indels -e 5 "
    f"--wildcard-file {expDir}{expName}_wcR2_L1_R2.txt "
    f"{ds_fq} "
    f"> {expDir}{expName}_L1_R2-ATAC.log"
)

dmuxL2 = (
    f"cutadapt "
    f"-g linker2={linker2} "
    f"-o \"{expDir}found-{expName}-"+"{name}"+f".fastq.gz\" "
    "--action=lowercase --cores 30 "
    "--no-indels -e 5 "
    f"--wildcard-file {expDir}{expName}_wcR2_L2_R2.txt "
    f"{ds_fq} "
    f"> {expDir}{expName}_L2_R2-ATAC.log"
)

logging.info(f"ds: {ds_cmd}\nl1: {dmuxL1}\nl2: {dmuxL2}")

# Run downsampling with seqtk, parse and filter barcodes with cutadapt ####
logging.info(f"Starting downsampling to {outReads} reads")
subprocess.run(ds_cmd, shell=True)
logging.info("Completed downsampling")

logging.info("Starting dmuxL1")
subprocess.run(dmuxL1, shell=True)
logging.info("Completed dmuxL1")

logging.info("Starting dmuxL2")
subprocess.run(dmuxL2, shell=True)
logging.info("Completed dmuxL2")

# For each dmux, process wildcard outputs.  Get read counts per mer, ####
# fractions, cumulative sums, map to plate info, identify expected ####
# barcode.  Add each table to countTableList, save to ####
# '{expDir}{eL}_{"-".join(bcPlate)}.csv'. ####
wildCardFileL1 = f"{expDir}{expName}_wcR2_L1_R2.txt"
wildCardFileL2 = f"{expDir}{expName}_wcR2_L2_R2.txt"
wcList = [wildCardFileL1, wildCardFileL2]
logging.info(wcList)

logFileL1 = f"{expDir}{expName}_L1_R2-ATAC.log"
logFileL2 = f"{expDir}{expName}_L2_R2-ATAC.log"
logList = [logFileL1, logFileL2]
logging.info(logList)

rowColList = ['spatialRow', 'spatialCol']
expList = [f"{expName}-{linkString}" for linkString in ["L1", "L2"]]
countTableList = []
maxToNinety = -1
outLines = []
for wc, bcPlate, logF, rc, eL in zip(
    wcList, plateList, logList, rowColList, expList
):
    logging.info(f"wcFile: {wc}\tplate: {bcPlate}\tlog: {logF}\tchannel: {rc}")
    countTable1s = pd.DataFrame()

    try:
        df = utils.load_wc_file(wc)
    except ValueError as e:
        logging.error(f"{e}")
        exit(0)

    uniqueCounts = df['8mer'].value_counts()
    countTable1 = pd.DataFrame(uniqueCounts)
    countTable1['fracCount'] = uniqueCounts/uniqueCounts.sum()

    countTable1s = countTable1.sort_values(
        by=['fracCount'], ascending=False
    ).copy()
    countTable1s['cumulative_sum'] = countTable1s['fracCount'].cumsum()
    numToNinety = (countTable1s['cumulative_sum'] <= 0.9).sum()
    try:
        pctFor50 = f"{countTable1s.iloc[:50,1].sum():.1%}"
        pctFor96 = f"{countTable1s.iloc[:96,1].sum():.1%}"
    except:
        pctFor50 = 'NaN'
        pctFor96 = 'NaN'

    # mark bcs that are expected from whitelist
    # intersection: what bcs are shared by both
    intSet = (
        set(countTable1s.index.tolist()) &
        set(merTable.loc[merTable.Plate.isin(bcPlate)].dropna(
            subset=['8mer'], inplace=False
        )['8mer'].tolist())
    )
    lenIntSet = len(intSet)

    # add column and set matches in set to True
    countTable1s['expectMer'] = False
    countTable1s['channel'] = -1
    countTable1s['merLabel'] = countTable1s.index

    for mer in list(intSet):
        countTable1s.loc[mer, 'expectMer'] = True

    countTableList.append(countTable1s)
    countTable1s.to_csv(f'{expDir}{eL}_{"-".join(bcPlate)}.csv', index=False)
    totalReadFromExpected = countTable1s['fracCount'][countTable1s['expectMer']].sum()
    results = subprocess.run(
        ['grep', 'Total reads\|Reads with', f'{logF}'],
        stdout=subprocess.PIPE
    ).stdout.decode('utf-8')
    tRead, aRead = utils.parseReadGrep(results)

    out = f'\n######### Info for {os.path.basename(wc)}\n'
    out = out + f'Total Reads: {tRead}  Reads with Adapter: {aRead}'
    out = out + f'\nThe number of unique strings in the 8mer column is {len(uniqueCounts)}.\nNinety percent (90%) of the reads come from total of {numToNinety} 8mers.'
    out = out + f"\nTotal of {lenIntSet} out of {len(merTable.loc[merTable.Plate.isin(bcPlate)].dropna(subset = ['8mer'], inplace=False))} expected 8-mers accounted for {totalReadFromExpected:.1%} of the reads"
    out = out + f"\nTop 50 8mers represent {pctFor50} fraction of reads\nTop 96 8mers represent {pctFor96} fraction of reads"

    outLines.append(out)
    logging.info(out)

    if numToNinety > maxToNinety:
        maxToNinety = numToNinety

# For each of cutadapt result (L1, L2, L1L2), identify hi/lo mers, save ####
# table as png if present, make barcode qc plot with hi/lo thresholds, ####
# and make barcode plate heatmap as png.
colAnalysis = 'fracCount'
picPaths = []
tablePicPaths = []
for bcp, ctb, eL in zip(plateList, countTableList, expList):
    expectedTable = ctb.copy()
    expectedTable['8mer'] = expectedTable.index
    expectedTable = expectedTable[expectedTable.channel != -1].sort_values(by=['channel'], ascending=True)
    expectedTable['hiWarn'] = False
    expectedTable['loWarn'] = False

    meanCount = expectedTable[colAnalysis].mean()
    upperCut = 2 * meanCount
    lowerCut = 0.5 * meanCount

    expectedTable.loc[expectedTable[colAnalysis] > upperCut, 'hiWarn'] = True
    expectedTable.loc[expectedTable[colAnalysis] < lowerCut, 'loWarn'] = True

    totalHiWarn = expectedTable['hiWarn'].sum()
    totalLoWarn = expectedTable['loWarn'].sum()
    totalMers = len(expectedTable[colAnalysis])
    out = f'\n######### Info for {eL}_{"-".join(bcp)}\n'
    out = out+f'Total hiWarn: {totalHiWarn}\tTotal loWarn: {totalLoWarn}'
    out = out+f'\nPct hiWarn: {(totalHiWarn/totalMers):.3f}\tPct loWarn: {(totalLoWarn/totalMers):.3f}'

    outLines.append(out)
    logging.info(out)

    # save hi/loWarn table
    if (totalHiWarn + totalLoWarn) > 0:  # Only export if there are hi or lows

        # Create a subset of the expectedTable where hiWarn and loWarn are true
        subset_expectedTable = expectedTable.loc[expectedTable['hiWarn'] | expectedTable['loWarn']]
        subset_expectedTable.to_csv(f'{expDir}{eL}_{"-".join(bcp)}_hiLoWarn.csv', index=False)

        # After generating the hi/loWarn table, save png of table
        hi_lo_warn_table_path = f'{expDir}{eL}_{"-".join(bcp)}_hiLoWarn.png'
        subset_expectedTable.loc[subset_expectedTable['hiWarn'] == True, 'hiOrLo'] = 'high'
        subset_expectedTable.loc[subset_expectedTable['loWarn'] == True, 'hiOrLo'] = 'low'
        expectColList = ['merLabel', 'plate', 'plateCol', 'plateRow', 'count', 'hiOrLo', 'fracCount']
        utils.save_table_as_image(
            subset_expectedTable[expectColList],
            hi_lo_warn_table_path
        )

        logging.info(f"Saved hi/loWarn table as image: {hi_lo_warn_table_path}")
        tablePicPaths.append(hi_lo_warn_table_path)

    # make bar plot (barcode QC plot)
    fig, ax = plt.subplots(figsize=(20, 5))
    expectedTable.plot(
        kind='bar', x='channel', y=colAnalysis, ax=ax, logy=False, legend=False
    )
    ax.set_xticks(expectedTable.channel)
    ax.set_xticklabels(expectedTable.merLabel, rotation=90)
    ax.axhline(y=meanCount, color='tab:gray', linestyle='--', linewidth=0.5, label=f'mean: {meanCount:.3f}')
    ax.axhline(y=meanCount/2, color='tab:red', linestyle='--', linewidth=0.5, label=f'half mean: {lowerCut:.3f}')
    ax.axhline(y=meanCount*2, color='tab:red', linestyle='--', linewidth=0.5, label=f'dbl mean: {upperCut:.3f}')
    fig.legend(loc="upper right")
    fig.tight_layout()

    # saving plot
    barplot_name = f'{expDir}{eL}_{"-".join(bcp)}_plot.png'
    logging.info(f"Saving plot {barplot_name}...")
    fig.savefig(barplot_name)
    picPaths.append(barplot_name)
    logging.info(f"Completed {barplot_name}")

    # Make barcode plate heatmaps
    for plate in bcp:
        tmpTable = expectedTable.loc[expectedTable.plate == plate].copy()
        pivoted = tmpTable.pivot(index="plateRow", columns="plateCol", values=colAnalysis)

        # Define the column and row labels
        cols = list(range(1, 13))
        rows = list("ABCDEFGH")
        pltMtx = np.zeros(shape=(len(rows), len(cols)))
        chMtx = np.zeros(shape=(len(rows), len(cols)))

        for i in range(len(rows)):
            for j in range(len(cols)):
                metVal = tmpTable.loc[(tmpTable['plateRow'] == rows[i]) & (tmpTable['plateCol'] == cols[j])]
                mval = metVal[colAnalysis].iloc[0] if len(metVal) != 0 else np.nan
                pltMtx[i, j] = mval
                chVal = tmpTable.loc[(tmpTable['plateRow'] == rows[i]) & (tmpTable['plateCol'] == cols[j])]
                cval = chVal['channel'].iloc[0] if len(chVal) != 0 else 0
                chMtx[i, j] = cval

        # Create a figure and an axes object
        fig, ax = plt.subplots()

        # Plot the heatmap using imshow
        im = ax.imshow(pltMtx, cmap="viridis", vmin=0)

        # Add a colorbar to show the scale
        cbar = ax.figure.colorbar(im, ax=ax, shrink=.75)

        # Set the tick labels and positions for the x and y axes
        ax.set_xticks(np.arange(len(cols)))
        ax.set_yticks(np.arange(len(rows)))
        ax.set_xticklabels(cols)
        ax.set_yticklabels(rows)

        # Rotate the x tick labels and set their alignment
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                 rotation_mode="anchor"
                 )

        # Loop over the data and create text annotations for each cell
        for i in range(len(rows)):
            for j in range(len(cols)):
                text = ax.text(j, i, int(chMtx[i, j]+1), ha="center",
                               va="center", color="w"
                               )

        # Set the title and adjust the layout
        ax.set_title(f'PlateMap: {plate}\n{eL}', fontsize=10)
        fig.tight_layout()

        # Save the plot
        platemap_name = f'{expDir}{eL}_{plate}_plate.png'
        logging.info(f"Saving {platemap_name}...")
        plt.savefig(platemap_name)
        picPaths.append(platemap_name)
        logging.info(f"Completed {platemap_name}")

# Make spatialTable.csv, spatial plot(s) ####
# if remote tissue position file was specified, use that, otherwise make one
if tissue_pos_path != "none":
    logging.info(f"tissue positions supplied: {tissue_pos_path}")
else:
    _, _, tissue_pos_path = utils.merList2BCList2TP(
        merPath, 'temp', outDirPath=expDir
    )
    logging.info(f"tissue positions generated: {tissue_pos_path}")

mtb = utils.wildcardSpatial(wildCardFileL1, wildCardFileL2, tissue_pos_path)
mtbPath = f"{expDir}{expName}_spatialTable.csv"
mtb.to_csv(mtbPath, index=False)

# Do some on/off tissue calcs ####
totalTix = len(mtb['onOff'])
totalOnTis = len(mtb.loc[mtb['onOff'] == 1]['onOff'])
totalOffTis = len(mtb.loc[mtb['onOff'] == 0]['onOff'])
countsOnTis = mtb.loc[mtb['onOff'] == 1]['count'].sum()
countsOffTis = mtb.loc[mtb['onOff'] == 0]['count'].sum()
countsPerTixOnTis = countsOnTis / totalOnTis

# generate density plot for on/off tissue pixels
density_plotPath = f'{expDir}{expName}_denseOnOff.png'

utils.create_density_plot(
    mtb,
    density_plotPath,
    'count',
    'onOff',
    log10=True,
    x_label='total counts',
    y_label='density'
)
picPaths.append(density_plotPath)

try:
    countsPerTixOffTis = countsOffTis / totalOffTis
except ZeroDivisionError:
    countsPerTixOffTis = 0

fractCPToffTiss = countsPerTixOffTis / countsPerTixOnTis
ratioOffvOn = countsOffTis / countsOnTis

logging.info(f"TixelStats: {totalTix=}, {totalOnTis=}, {totalOffTis=}, \
             {countsOnTis=}, {countsOffTis=}, {ratioOffvOn=}, \
             {countsPerTixOnTis=}, {countsPerTixOffTis=}, {fractCPToffTiss=}")

valList = [totalTix, totalOnTis, totalOffTis, countsOnTis,
           countsOffTis, countsPerTixOnTis, countsPerTixOffTis,
           fractCPToffTiss, ratioOffvOn]

nameList = ['totalTix', 'totalOnTis', 'totalOffTis', 'countsOnTis',
            'countsOffTis', 'countsPerTixOnTis', 'countsPerTixOffTis',
            'fractCPToffTiss', 'ratioOffvOn']

mTable = utils.variables_to_dataframe(valList, nameList)
mTablePath = f"{expDir}{expName}_onOffTissueTable_mqc.csv"
mTable.to_csv(mTablePath, index=False)

# write html with plate/plot figures, append table pics to the end of ####
#  picPaths ####
picPaths = picPaths + tablePicPaths
reportNoteHTML = f"<br><b>Reads with L1:</b><br>{tRead}<br><b>Reads with L1&L2:</b><br>{aRead}"
utils.generate_html_with_embedded_images(
    picPaths, expDir, expName, noteHTML=reportNoteHTML
)

# Make pareto chart of barcode abundances, save as _output.pdf ####
with PdfPages(f'{expDir}{expName}_output.pdf') as pdf:
    for tb, wc, ol in zip(countTableList, wcList, outLines):
        # define subplots
        startI = 0
        endI = maxToNinety + int(maxToNinety / 2)
        fig, ax = plt.subplots(figsize=(20, 10))

        ax.plot(tb['fracCount'][startI:endI], marker='x',
                linewidth=1, label='Fraction total 8mer')

        # define second y-axis that shares x-axis with current plot
        ax2 = ax.twinx()

        # add second line to plot
        ax2.plot(
            tb['cumulative_sum'][startI:endI], marker='.', linewidth=1,
            color='tab:orange', label='Cumlative Fraction'
        )

        # add second y-axis label
        labels0 = list(tb.index[startI:endI])
        labels2 = list(tb.index[startI:endI] + "_" + tb.channel[startI:endI].astype(str).str.zfill(2))

        chanLabels = list("_" + (tb.channel[startI:endI] + 1).astype(str).str.zfill(2))
        labels3 = [f"{lb}{cb}" for lb, cb in zip(labels0, chanLabels)]

        labels = labels3
        ax.set_xticks(labels0)
        ax.set_xticklabels(labels, rotation=90)

        s = tb.expectMer[startI:endI]
        makeGreen = np.where(s)[0].tolist()
        for mg in makeGreen:
            ax.get_xticklabels()[mg].set_color("green")

        ax2.axhline(y=.9, color='tab:gray', linestyle='--', linewidth=0.5, label='90 percent')
        fig.legend(loc="center right")
        plt.title(f"plot: {wc}")
        plt.annotate(ol, xy=(0.6, 0.6), xycoords='axes fraction', ha='left', va='center', ma='left')
        fig.subplots_adjust(bottom=0.2)
        plt.close()

        pdf.savefig(fig)
