from __future__ import annotations

import logging
import pandas as pd
import re

from importlib.resources import files as resource_files
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


def contains_acgt_word(input_list):
    '''Function to check for 8-character word made up of A, C, G, T in a list
    and return indices.'''
    pattern = re.compile(r'^[ACGTN]{8,}$')
    return [
        index for index, item in enumerate(input_list) if pattern.search(item)
    ]


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
