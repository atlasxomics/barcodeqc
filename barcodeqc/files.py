import logging
import pandas as pd

from pathlib import Path

import barcodeqc.utils as utils

logger = logging.getLogger(__name__)


class BarcodeFileError(ValueError):
    """Raised when a positions or barcode file fails validation."""


class WildcardFileError(ValueError):
    """Raised when a positions or barcode file fails validation."""


def load_wc_file(wcPath):
    '''Loads a space-delimited text file containing parsed read names from a
    fastq file.  Checks to see if one of the rows contains an 8mer barcode.  If
    it does, assumes the next column contains the read and and returns a
    DataFrame, with the columns renamed.
    '''
    df = pd.read_csv(wcPath, header=None, sep=' ')
    colNames = list(range(0, len(df.columns)))
    testLine = list(df.iloc[5, :])

    if len(merCol := utils.contains_acgt_word(testLine)) != 1:
        raise WildcardFileError(f"No 8-mer column matches found: {testLine=}")

    colNames[merCol[0]] = '8mer'
    colNames[merCol[0] + 1] = 'readName'
    df.columns = colNames

    return df


def open_barcode_file(bc_path: Path) -> pd.DataFrame:
    """Open a barcode file and ensure proper format.
    Expected tsv of (three columns):
      0: sequence (string, contains an 16-mer of A/T/C/G/N)
      2: row position (non-negative integer)
      3: column position (non-negative integer)

    Parameters
    -------
    bc_path : Path
        Path to barcode file.

    Returns
    -------
    pd.DataFrame
        Validated positions table with columns:
        ['sequence', 'row', 'col']
    """

    barcodes = pd.read_csv(bc_path, header=0)

    if barcodes.shape[1] < 3:
        raise BarcodeFileError(
            f"Positions file must have at least 3 columns; "
            f"found {barcodes.shape[1]}"
        )

    barcodes = barcodes.iloc[:, :3]

    # must be string-like
    seq = barcodes["sequence"]
    if not seq.map(lambda x: isinstance(x, str)).all():
        bad = barcodes.loc[~seq.map(lambda x: isinstance(x, str)), "barcodes"]
        raise BarcodeFileError(
            f"Non-string barcodes detected (example: {bad.iloc[0]!r})"
            "Please ensure you are using the correct tissue_postions file"
        )

    # must be exactly 8 chars of A,T,C,G,N
    pattern = r"^[ATCGN]{8}$"
    invalid = ~barcodes["sequence"].astype(str).str.contains(pattern)

    if invalid.any():
        bad = seq[invalid].unique()[:5]
        raise BarcodeFileError(
            f"Invalid barcodes detected (expected 8-mer of A/T/C/G/N). "
            f"Examples: {bad}"
            "Please ensure you are using the correct tissue_postions file"
        )

    for c in ["row", "col"]:
        vals = pd.to_numeric(barcodes[c], errors="coerce")

        # non-numeric values
        if vals.isna().any():
            bad = barcodes.loc[vals.isna(), c].unique()[:5]
            raise BarcodeFileError(
                f"{c} column contains non-numeric values. Examples: {bad}"
            )

        # must be integer-valued
        if (vals % 1 != 0).any():
            bad = vals[vals % 1 != 0].unique()[:5]
            raise BarcodeFileError(
                f"{c} column must contain integer values. Examples: {bad}"
            )

        # must be non-negative
        if (vals < 0).any():
            bad = vals[vals < 0].unique()[:5]
            raise BarcodeFileError(
                f"{c} column contains negative values. Examples: {bad}"
            )

        # assign cleaned column back
        barcodes[c] = vals.astype(int)

    return barcodes


def open_positions_file(position_path: Path) -> pd.DataFrame:
    """Open a tissue_positions_file and validate barcode fidelity.
    Expected format (first 4 columns):
      0: barcode (string, contains an 16-mer of A/T/C/G/N; optional
        trailing '-1')
      1: on_off flag (0 or 1)
      2: row position (non-negative integer)
      3: column position (non-negative integer)

    Returns
    -------
    pd.DataFrame
        Validated positions table with columns:
        ['barcodes', 'on_off', 'row', 'col']
    """

    positions = pd.read_csv(position_path, header=None)

    if positions.shape[1] < 4:
        raise BarcodeFileError(
            f"Positions file must have at least 4 columns; "
            f"found {positions.shape[1]}"
            "Please ensure you are using the correct tissue_postions file"
        )

    positions = positions.iloc[:, :4]
    positions.columns = ["barcodes", "on_off", "row", "col"]

    # must be string-like
    bc = positions["barcodes"]
    if not bc.map(lambda x: isinstance(x, str)).all():
        bad = positions.loc[~bc.map(lambda x: isinstance(x, str)), "barcodes"]
        raise BarcodeFileError(
            f"Non-string barcodes detected (example: {bad.iloc[0]!r})"
            "Please ensure you are using the correct tissue_postions file"
        )

    # remove trailing "-1" only if it is at the end
    positions["barcodes"] = (
        positions["barcodes"]
        .astype(str)
        .str.replace(r"-1$", "", regex=True)
    )

    # must be exactly 8 chars of A,T,C,G,N
    pattern = r"^[ATCGN]{16}$"
    invalid = ~positions["barcodes"].astype(str).str.contains(pattern)

    if invalid.any():
        bad = bc[invalid].unique()[:5]
        raise BarcodeFileError(
            f"Invalid barcodes detected (expected 8-mer of A/T/C/G/N). "
            f"Examples: {bad}"
            "Please ensure you are using the correct tissue_postions file"
        )

    # must be exactly 0 or 1
    on_off = pd.to_numeric(positions["on_off"], errors="coerce")
    invalid = ~on_off.isin([0, 1])
    if invalid.any():
        bad = on_off[invalid].unique()[:5]
        raise BarcodeFileError(
            f"on_off column must contain only 0 or 1. Found: {bad}"
            "Please ensure you are using the correct tissue_postions file"
        )

    # assign cleaned column back
    positions["on_off"] = on_off.astype(int)

    for c in ["row", "col"]:
        vals = pd.to_numeric(positions[c], errors="coerce")

        # non-numeric values
        if vals.isna().any():
            bad = positions.loc[vals.isna(), c].unique()[:5]
            raise BarcodeFileError(
                f"{c} column contains non-numeric values. Examples: {bad}"
                "Please ensure you are using the correct tissue_postions file"
            )

        # must be integer-valued
        if (vals % 1 != 0).any():
            bad = vals[vals % 1 != 0].unique()[:5]
            raise BarcodeFileError(
                f"{c} column must contain integer values. Examples: {bad}"
                "Please ensure you are using the correct tissue_postions file"
            )

        # must be non-negative
        if (vals < 0).any():
            bad = vals[vals < 0].unique()[:5]
            raise BarcodeFileError(
                f"{c} column contains negative values. Examples: {bad}"
                "Please ensure you are using the correct tissue_postions file"
            )

        # assign cleaned column back
        positions[c] = vals.astype(int)

    return positions
