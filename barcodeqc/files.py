import pandas as pd

from pathlib import Path


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
        raise ValueError(
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
        raise ValueError(
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
        raise ValueError(
            f"Invalid barcodes detected (expected 8-mer of A/T/C/G/N). "
            f"Examples: {bad}"
            "Please ensure you are using the correct tissue_postions file"
        )

    # must be exactly 0 or 1
    on_off = pd.to_numeric(positions["on_off"], errors="coerce")
    invalid = ~on_off.isin([0, 1])
    if invalid.any():
        bad = on_off[invalid].unique()[:5]
        raise ValueError(
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
            raise ValueError(
                f"{c} column contains non-numeric values. Examples: {bad}"
                "Please ensure you are using the correct tissue_postions file"
            )

        # must be integer-valued
        if (vals % 1 != 0).any():
            bad = vals[vals % 1 != 0].unique()[:5]
            raise ValueError(
                f"{c} column must contain integer values. Examples: {bad}"
                "Please ensure you are using the correct tissue_postions file"
            )

        # must be non-negative
        if (vals < 0).any():
            bad = vals[vals < 0].unique()[:5]
            raise ValueError(
                f"{c} column contains negative values. Examples: {bad}"
                "Please ensure you are using the correct tissue_postions file"
            )

        # assign cleaned column back
        positions[c] = vals.astype(int)

    return positions
