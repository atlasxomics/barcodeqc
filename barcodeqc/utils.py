import base64
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
import seaborn as sns
import warnings

from jinja2 import Template
from matplotlib.ticker import FuncFormatter
from pathlib import Path

from barcodeqc.files import open_positions_file

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s",
    level=logging.INFO
)

repo_root = os.path.abspath(os.path.dirname(__file__))

BARCODE_PATHS = {
    "bc50": {
        "whitelist": f"{repo_root}/whitelists/bc50.txt",
        "positions": f"{repo_root}/position_files/x50_all_tissue_positions_list.csv"
    },
    "bc96": {
        "whitelist": f"{repo_root}/whitelists/bc96.txt",
        "positions": f"{repo_root}/position_files/x96_all_tissue_positions_list.csv"
    },
    "fg96": {
        "whitelist": f"{repo_root}/whitelists/bc96FG_11DEC.txt",
        "positions": f"{repo_root}/position_files/xfg96_11DEC_alltissue_positions_list.csv"
    },
    "bc220": {
        "whitelist": f"{repo_root}/whitelists/bc220_25APR.txt",
        "positions": f"{repo_root}/position_files/xbc220_25APR_alltissue_positions_list.csv"
    },
    "bc220_05-OCT": {
        "whitelist": f"{repo_root}/whitelists/bc220_05OCT.txt",
        "positions": f"{repo_root}/position_files/xbc220_05OCT_alltissue_positions_list.csv"
    },
    "bc220_20-MAY": {
        "whitelist": f"{repo_root}/whitelists/bc220_20MAY.txt",
        "positions": f"{repo_root}/position_files/xbc220-20MAY_alltissue_positions_list.csv"
    }
}


def concatenate_2pairs(list1, list2):
    '''Used to combine barcode A 8mer and barcode B 8mer'''
    return [b + a for a in list1 for b in list2]


def concatenate_2pairsChan(mer1, mer2, chan1, chan2):
    '''Combined barcode 8mers and row/column indices'''
    merPairs = [b + a for a in mer1 for b in mer2]
    chanPairs = [[int(c1), int(c2)] for c1 in chan1 for c2 in chan2]

    return merPairs, chanPairs


def contains_acgt_word(input_list):
    '''Function to check for 8-character word made up of A, C, G, T in a list
    and return indices.'''
    pattern = re.compile(r'^[ACGTN]{8,}$')
    return [
        index for index, item in enumerate(input_list) if pattern.search(item)
    ]


def create_density_plot(
    dataframe, outPath, data_column, group_column,
    log10=True, x_label='', y_label='', plotTitle=''
):
    """Generate and save a density plot; used downstream to plot the
    distribution of counts per barcode (_denseOnOff.png').
    """
    plt.figure(figsize=(6, 6))

    # Use the original palette for up to two groups
    palette_full = sns.color_palette("bright", 10)
    palette = [palette_full[0], palette_full[3]]

    # Get unique, non-null groups
    groups = [g for g in dataframe[group_column].dropna().unique()
              if len(dataframe[dataframe[group_column] == g]) > 0]

    if len(groups) == 0:
        logging.warning("No non-empty groups to plot in density plot.")
        return

    # Plot each non-empty group with the correct color
    for idx, group in enumerate(groups):
        group_df = dataframe[dataframe[group_column] == group]
        sns.kdeplot(
            data=group_df,
            x=data_column,
            bw_adjust=1.5,
            log_scale=log10,
            fill=True,
            color=palette[idx],
            label=f"{group}"
        )

    if log10:
        x_label = x_label + " (log scaled)"

    plt.title(plotTitle, fontsize=16, fontweight='bold')
    plt.xlabel(x_label, fontsize=14, fontweight='bold')
    plt.ylabel(y_label, fontsize=14, fontweight='bold')
    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:.2f}'))
    plt.legend(loc='upper right', frameon=False)
    plt.xlim(0, 10000)
    plt.savefig(f'{outPath}', dpi=300)
    plt.close()


def create_heatmap(
    df,
    outPath,
    log10=True,
    axesOff=False,
    countCol='nCellHits',
    vmin=None,
    vmax=None,
    colorMap='RdYlBu_r'
):
    '''Create spatial heatmap arranged by row/column, colored by read counts'''

    # Pivot the dataframe to create a matrix
    matrix = df.pivot(index='row', columns='col', values=countCol)
    cbLbl = "counts"

    # Take the log10 of the values in the matrix
    with warnings.catch_warnings():
        # if there is a zero, log10 throws a warning.
        # ignoring that warning
        warnings.simplefilter("ignore", RuntimeWarning)
        if log10:
            matrix = np.log10(matrix)
            cbLbl = "log10 of counts"
            if vmin is not None:
                vmin = np.log10(vmin)
            if vmax is not None:
                vmax = np.log10(vmax)

    # Create the heatmap
    fig, ax = plt.subplots(figsize=(15, 12))

    im = ax.imshow(matrix, vmin=vmin, vmax=vmax, cmap=colorMap)

    # Add labels
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_yticks(range(len(matrix.index)))
    ax.set_xticklabels(matrix.columns, rotation=90, fontsize=6)
    ax.set_yticklabels(matrix.index, fontsize=6)

    if axesOff:
        plt.axis('off')

    # Add colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    # cbar.set_label(cbLbl, rotation=270)  # changed to add labelpad
    cbar.set_label(cbLbl, rotation=270, labelpad=15)

    plt.tight_layout()  # Automatically adjust subplot parameters

    # save the plot
    # plt.savefig(outPath)  # change to reduce whitespace
    plt.savefig(outPath, bbox_inches='tight', pad_inches=0.1)

    logging.info(f"create_heatmap: Saving to {outPath}")
    plt.close()

    return matrix


def generate_html_with_embedded_images(
    picList, outPath, sampleName, fileTag='bcQC', noteHTML=''
):
    # Template for our HTML document

    html_template = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>AtlasXomics Barcode QC Report</title>
            <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css">
            <style>
                .fixed-sidebar {
                    position: fixed;
                    #top: 60;
                    height: 100%;
                    overflow-y: auto;
                }
                body {
                    #padding-top: 70px;
                    #top: 70px;
                }

            </style>


        </head>
        <body>
        <!---
        <header>
            <nav class="navbar navbar-dark bg-dark fixed-top">
                <a class="navbar-brand" href="#">SampleName: {{ sampleName }}</a>
            </nav>
        </header>
        --!>
            <div class="container-fluid">
                <div class="row">
                    <div class="col-3 fixed-sidebar">
                        <div class="list-group">
                        <a href="#" class="list-group-item list-group-item-action"><b>SampleName:</b><br>{{ sampleName }}{{ noteHTML }}</a>
                            {% for file, data_uri in files %}
                            <a href="#{{ file }}" class="list-group-item list-group-item-action">{{ file }}</a>
                            {% endfor %}
                        </div>
                    </div>
                    <div class="col-9 offset-3">
                        {% for file, data_uri in files %}
                        <h3 id="{{ file }}">{{ file }}</h3>
                        <img src="{{ data_uri }}" class="img-fluid" alt="{{ file }}">
                        {% endfor %}
                    </div>
                </div>
            </div>
            <script src="https://code.jquery.com/jquery-3.3.1.slim.min.js"></script>
            <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/js/bootstrap.min.js"></script>
        </body>
        </html>
        """

    # Function to get base64 encoded string for an image file
    def get_image_data_uri(filename):
        with open(filename, 'rb') as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
            return f'data:image/png;base64,{encoded_string}'

    # Get all .png files in the directory and their base64 encoded strings
    files_data_uri = [
        (os.path.basename(file), get_image_data_uri(file)) for file in picList
    ]

    # Render the HTML with the embedded images
    template = Template(html_template)
    html_content = template.render(
        files=files_data_uri, sampleName=sampleName, noteHTML=noteHTML
    )

    # Write the HTML content to a file
    with open(os.path.join(outPath, f'{sampleName}_{fileTag}_report.html'), 'w') as f:
        f.write(html_content)


def load_wc_file(wcPath):
    '''Loads a space-delimited text file containing parsed read names from a
    fastq file.  Checks to see if one of the rows contains an 8mer barcode.  If
    it does, assumes the next column contains the read and and returns a
    DataFrame, with the columns renamed.'''
    df1 = pd.read_csv(wcPath, header=None, sep=' ')
    colNames = list(range(0, len(df1.columns)))
    testLine = list(df1.iloc[5, :])

    if len(merCol := contains_acgt_word(testLine)) != 1:
        raise ValueError(f"No 8-mer column matches found: {testLine=}")

    colNames[merCol[0]] = '8mer'
    colNames[merCol[0] + 1] = 'readName'
    df1.columns = colNames
    return df1


def parse_read_log(log_path: str) -> tuple[str, str]:
    text = Path(log_path).read_text()

    tot = re.search(r"Total reads:\s*(\d+)", text)
    adapt = re.search(r"Reads with.*?:\s*(\d+)", text)

    if not tot or not adapt:
        raise ValueError("Expected read counts not found")

    return tot.group(1), adapt.group(1)


def save_table_as_image(df, file_path):
    fig, ax = plt.subplots(figsize=(12, 8))  # set size frame
    ax.axis('tight')
    ax.axis('off')
    table = ax.table(
        cellText=df.values,
        colLabels=df.columns,
        cellLoc='center',
        loc='center'
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.2)  # may need to adjust depending on your data
    plt.savefig(file_path, bbox_inches='tight')
    plt.close()


def variables_to_dataframe(vlist, nlist):
    """
    Converts a list of variables to a pandas DataFrame.
    Used to make onOffTissueTable_mqc.csv.

    Args:
        variables (list): List of variables (can be any data type).

    Returns:
        pd.DataFrame: DataFrame with two columns: 'Variable' and 'Value'.
    """
    data = []
    for var_name, var_value in zip(nlist, vlist):
        data.append({'Metric': var_name, 'Value': var_value})

    df = pd.DataFrame(data)
    return df


def wildcardSpatial(wcL1File, wcL2File, tissuePosnFile):
    '''Merges wild-card output files from cutadapt, then merges with tissue
    position file to create base for _spatialTable.csv.  Returns Dataframe.
    '''
    # read read2 L1 wc list
    try:
        wcL1tbl = load_wc_file(wcL1File)
    except ValueError as e:
        logging.error(f"{e}")
        exit(0)
    wcL1tbl['8mer_L1'] = wcL1tbl['8mer']

    # read read2 L2 wc List
    try:
        wcL2tbl = load_wc_file(wcL2File)
    except ValueError as e:
        logging.error(f"{e}")
        exit(0)
    wcL2tbl['8mer_L2'] = wcL2tbl['8mer']

    # merge on second column
    mergedTable8Mers = pd.merge(wcL1tbl, wcL2tbl, on='readName')

    # concat 8-mers. In the read, B is first (associated with L2)
    mergedTable8Mers['16mer'] = mergedTable8Mers['8mer_L2'] + mergedTable8Mers['8mer_L1']

    # count table of 16mers
    merCount = mergedTable8Mers['16mer'].value_counts()
    countTable16mer = pd.DataFrame(merCount)

    # Read in tissue position file.
    tissue_positions = open_positions_file(tissuePosnFile, header=None)

    # merge tissue position and countTable16mer
    mergedTablePosition = pd.merge(
        countTable16mer,
        tissue_positions,
        left_on='16mer',
        right_on='barcodes',
        how='outer'
    )

    # remove all 16mers that are NOT in the tissue position file
    mergedTablePosition.dropna(subset=['row', 'col', 'on_off'], inplace=True)

    return mergedTablePosition
