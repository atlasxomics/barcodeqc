from pathlib import Path

from importlib.resources import files as files


PACKAGE_DIR = Path(__file__).parent
DATA_DIR = files("barcodeqc") / "data"

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
