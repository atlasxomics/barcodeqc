from __future__ import annotations

import base64
import logging
import os
import pandas as pd
import re
import sys
from datetime import datetime

from importlib.resources import files as resource_files
from jinja2 import Template
from pathlib import Path

import barcodeqc.files as files

logger = logging.getLogger(__name__)

PACKAGE_DIR = Path(__file__).parent
DATA_DIR = resource_files("barcodeqc") / "data"

DATA_DIRS = {
    "positions": DATA_DIR / "position_files",
    "barcodes": DATA_DIR / "barcode_files",
}

BARCODE_PATHS = {
    "bc50": {
        "positions": DATA_DIRS["positions"] / "x50_all_tissue_positions_list.csv",
        "bca": DATA_DIRS["barcodes"] / "barcode_A" / "merList50.tsv",
        "bcb": DATA_DIRS["barcodes"] / "barcode_B" / "merList50.tsv"
    },
    "bc96": {
        "positions": DATA_DIRS["positions"] / "x96_all_tissue_positions_list.csv",
        "bca": DATA_DIRS["barcodes"] / "barcode_A" / "merList96.tsv",
        "bcb": DATA_DIRS["barcodes"] / "barcode_B" / "merList96.tsv"
    },
    "fg96": {
        "positions": DATA_DIRS["positions"] / "xfg96_11DEC_alltissue_positions_list.csv",
        "bca": DATA_DIRS["barcodes"] / "barcode_A" / "merListfg96.tsv",
        "bcb": DATA_DIRS["barcodes"] / "barcode_B" / "merListfg96.tsv"
    },
    "bc220": {
        "positions": DATA_DIRS["positions"] / "xbc220_25APR_alltissue_positions_list.csv",
        "bca": DATA_DIRS["barcodes"] / "barcode_A" / "merList220_25-APR.tsv",
        "bcb": DATA_DIRS["barcodes"] / "barcode_B" / "merList220_25-APR.tsv"
    },
    "bc220_05-OCT": {
        "positions": DATA_DIRS["positions"] / "xbc220_05OCT_alltissue_positions_list.csv",
        "bca": DATA_DIRS["barcodes"] / "barcode_A" / "merList220_05-OCT.tsv",
        "bcb": DATA_DIRS["barcodes"] / "barcode_B" / "merList220_05-OCT.tsv"
    },
    "bc220_20-MAY": {
        "positions": DATA_DIRS["positions"] / "xbc220-20MAY_alltissue_positions_list.csv",
        "bca": DATA_DIRS["barcodes"] / "barcode_A" / "merList220_20-MAY.tsv",
        "bcb": DATA_DIRS["barcodes"] / "barcode_B" / "merList220_20-MAY.tsv"
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


def make_spatial_table(wcL1File, wcL2File, tissuePosnFile):
    '''Merges wild-card output files from cutadapt, then merges with tissue
    position file to create base for _spatialTable.csv.  Returns Dataframe.
    '''
    # read read2 L1 wc list
    try:
        wcL1tbl = files.load_wc_file(wcL1File)
    except ValueError as e:
        logger.error(f"{e}")
        exit(0)
    wcL1tbl['8mer_L1'] = wcL1tbl['8mer']

    # read read2 L2 wc List
    try:
        wcL2tbl = files.load_wc_file(wcL2File)
    except ValueError as e:
        logger.error(f"{e}")
        exit(0)
    wcL2tbl['8mer_L2'] = wcL2tbl['8mer']

    # merge on second column
    mergedTable8Mers = pd.merge(wcL1tbl, wcL2tbl, on='readName')

    # concat 8-mers. In the read, B is first (associated with L2)
    mergedTable8Mers['16mer'] = mergedTable8Mers['8mer_L2'] + mergedTable8Mers['8mer_L1']

    # count table of 16mers
    merCount = mergedTable8Mers['16mer'].value_counts()
    countTable16mer = pd.DataFrame(merCount)
    countTable16mer = countTable16mer.reset_index()
    countTable16mer.columns = ['16mer', 'count']
    countTable16mer['16mer'] = countTable16mer['16mer'].astype(str)

    # Read in tissue position file.
    tissue_positions = files.open_positions_file(tissuePosnFile)
    tissue_positions['barcodes'] = tissue_positions['barcodes'].astype(str)

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


def parse_read_log(log_path: str) -> tuple[str, str]:
    text = Path(log_path).read_text()

    tot = re.search(r"Total reads processed:\s*([\d,]+)", text)
    adapt = re.search(r"Reads with adapters:\s*([\d,]+)", text)

    if not tot or not adapt:
        raise ValueError("Expected read counts not found")

    total_reads = tot.group(1).replace(",", "")
    adapter_reads = adapt.group(1).replace(",", "")
    return total_reads, adapter_reads


def setup_logging(
    *,
    log_file: str | Path | None = None,
    stdout_level: int = logging.INFO,
    file_level: int = logging.DEBUG,
    log_dir: Path | None = None,
) -> None:
    """
    Configure application-wide logging.

    This should be called exactly once, at CLI startup.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    formatter = logging.Formatter(
        "%(levelname)s - %(asctime)s - %(message)s",
        "%Y-%m-%d %H:%M:%S",
    )

    # stdout
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(stdout_level)
    sh.setFormatter(formatter)
    logger.addHandler(sh)

    # optional file
    if log_file is not None:
        if log_dir is not None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            log_file = Path(log_dir) / f"{Path(log_file).stem}_{timestamp}.log"
        else:
            log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)

        fh = logging.FileHandler(log_file)
        fh.setLevel(file_level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
