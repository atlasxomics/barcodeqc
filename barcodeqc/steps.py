from __future__ import annotations

import logging
import subprocess

from pathlib import Path

import pandas as pd

from barcodeqc.config import QCConfig
import barcodeqc.files as files

logger = logging.getLogger(__name__)


def build_spatial_table(
    wc_linker1: Path,
    wc_linker2: Path,
    tissue_position_file: Path,
    output_dir: Path,
) -> tuple[pd.DataFrame, Path]:
    spatial_table = make_spatial_table(
        wc_linker1, wc_linker2, tissue_position_file
    )
    spatial_table_path = output_dir / "spatialTable.csv"
    spatial_table.to_csv(spatial_table_path, index=False)
    return spatial_table, spatial_table_path


def build_count_table(
    wc_path: Path,
    bcl: pd.DataFrame,
    row_col: str,
) -> tuple[pd.DataFrame, pd.Series, set[str], int, str, str, pd.Series]:

    whitelist = bcl["sequence"]
    wc_df = files.load_wc_file(wc_path)

    unique_counts = wc_df["8mer"].value_counts()
    count_table = unique_counts.to_frame(name="count")
    count_table["frac_count"] = unique_counts / unique_counts.sum()

    count_table = count_table.sort_values(by=["frac_count"], ascending=False)
    count_table["cumulative_sum"] = count_table["frac_count"].cumsum()
    num_to_ninety = (count_table["cumulative_sum"] <= 0.9).sum()

    try:
        pct_for_50 = f"{count_table.iloc[:50, 1].sum():.1%}"
        pct_for_96 = f"{count_table.iloc[:96, 1].sum():.1%}"
    except Exception:
        pct_for_50 = "NaN"
        pct_for_96 = "NaN"

    expected_bcs = set(count_table.index.tolist()) & set(whitelist)

    count_table = count_table.join(
        bcl.set_index("sequence")[["row", "col"]],
        how="left",
    )
    count_table["sequence"] = count_table.index

    count_table["channel"] = count_table[row_col]
    count_table["expectMer"] = count_table["sequence"].isin(expected_bcs)

    return (
        count_table,
        unique_counts,
        expected_bcs,
        num_to_ninety,
        pct_for_50,
        pct_for_96,
        whitelist
    )


def compute_hi_lo_qc(
    count_table: pd.DataFrame,
    channel_col: str = "channel",
    frac_col: str = "frac_count",
    expect_col: str = "expectMer",
) -> tuple[pd.DataFrame, int, int, int]:
    bc_table = count_table[count_table[expect_col]].sort_values(
        by=[channel_col], ascending=True
    )

    bc_table = bc_table.copy()
    bc_table["hiWarn"] = False
    bc_table["loWarn"] = False

    mean = bc_table[frac_col].mean()
    upper_cut = 2 * mean
    lower_cut = 0.5 * mean

    bc_table.loc[bc_table[frac_col] > upper_cut, "hiWarn"] = True
    bc_table.loc[bc_table[frac_col] < lower_cut, "loWarn"] = True

    total_hi_warn = int(bc_table["hiWarn"].sum())
    total_lo_warn = int(bc_table["loWarn"].sum())
    total_mers = int(len(bc_table))

    return bc_table, total_hi_warn, total_lo_warn, total_mers


def compute_onoff_metrics(spatial_table: pd.DataFrame) -> pd.DataFrame:
    on_df = spatial_table.loc[spatial_table["on_off"] == 1]
    off_df = spatial_table.loc[spatial_table["on_off"] == 0]

    total_pix = len(spatial_table)
    total_on = len(on_df)
    total_off = len(off_df)
    counts_on = on_df["count"].sum()
    counts_off = off_df["count"].sum()
    counts_per_pix_on = counts_on / total_on if total_on > 0 else 0
    counts_per_pix_off = counts_off / total_off if total_off > 0 else 0
    frac_per_pix_off_on = (
        counts_per_pix_off / counts_per_pix_on
        if counts_per_pix_on > 0
        else 0
    )
    ratio_off_on = counts_off / counts_on if counts_on > 0 else 0

    return pd.DataFrame(
        {
            "metric": [
                "total_pix",
                "total_on",
                "total_off",
                "counts_on",
                "counts_off",
                "ratio_off_on",
                "counts_per_pix_on",
                "counts_per_pix_off",
                "frac_per_pix_off_on",
            ],
            "value": [
                total_pix,
                total_on,
                total_off,
                counts_on,
                counts_off,
                ratio_off_on,
                counts_per_pix_on,
                counts_per_pix_off,
                frac_per_pix_off_on,
            ],
        }
    )


def ensure_output_dir(output_dir: Path) -> None:
    if not output_dir.exists():
        logger.debug(f"Output dir {output_dir} not detected, making...")
        output_dir.mkdir(parents=True, exist_ok=True)


def linker_conservation_status(
    total_reads: int,
    adapter_reads: int,
    pass_threshold: float = 0.7,
) -> tuple[str, float]:
    if total_reads <= 0:
        return "CAUTION", 0.0
    pct = adapter_reads / total_reads
    status = "PASS" if pct >= pass_threshold else "CAUTION"
    return status, pct


def barcode_check_status(
    count_table: pd.DataFrame,
    expect_col: str = "expectMer",
    frac_col: str = "frac_count",
    top_n: int = 100,
) -> tuple[str, int]:
    top = count_table.sort_values(by=[frac_col], ascending=False).head(top_n)
    unexpected = int((~top[expect_col]).sum())
    status = "PASS" if unexpected == 0 else "CAUTION"
    return status, unexpected


def lane_status(
    bc_table: pd.DataFrame,
    flag_col: str,
    channel_col: str = "channel",
) -> str:
    flagged = (
        bc_table.loc[bc_table[flag_col], channel_col]
        .dropna()
        .astype(int)
        .unique()
        .tolist()
    )
    if not flagged:
        return "PASS"
    flagged_sorted = sorted(flagged)
    adjacent = any(
        b - a == 1 for a, b in zip(flagged_sorted, flagged_sorted[1:])
    )
    return "CONTACT SUPPORT" if adjacent else "ACTION REQUIRED"


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


def run_subsample(config: QCConfig, output_dir: Path) -> Path:
    ds_path = output_dir / f"ds_{config.sample_reads}.fastq.gz"
    seqtk_cmd = [
        "seqtk",
        "sample",
        "-s",
        str(config.random_seed),
        str(config.r2_path),
        str(config.sample_reads),
    ]
    gzip_cmd = ["gzip"]
    logger.info("Running subsample ")
    logger.debug("seqtk cmd: %s", seqtk_cmd)
    logger.debug("gzip cmd: %s > %s", gzip_cmd, ds_path)
    with open(ds_path, "wb") as out_f:
        seqtk_proc = subprocess.Popen(seqtk_cmd, stdout=subprocess.PIPE)
        gzip_proc = subprocess.run(
            gzip_cmd,
            stdin=seqtk_proc.stdout,
            stdout=out_f,
            check=True,
        )
        if seqtk_proc.stdout is not None:
            seqtk_proc.stdout.close()
        seqtk_return = seqtk_proc.wait()
        if seqtk_return != 0:
            raise subprocess.CalledProcessError(seqtk_return, seqtk_cmd)
    logger.info("Completed subsampling")
    return ds_path


def run_cutadapt(
    ds_path: Path,
    output_dir: Path,
    cores: int = 30,
) -> tuple[Path, Path, Path, Path]:
    linker1 = "NNNNNNNNGTGGCCGATGTTTCGCATCGGCGTACGACT"
    linker2 = "NNNNNNNNATCCACGTGCTTGAGAGGCCAGAGCATTCG"

    wc_linker1 = output_dir / "cutadapt_wc_L1.txt"
    log_linker1 = output_dir / "cutadapt_L1.log"
    dmuxL1 = [
        "cutadapt",
        "-g",
        f"linker1={linker1}",
        "-o",
        "/dev/null",
        "--action=lowercase",
        "--cores",
        str(cores),
        "--no-indels",
        "-e",
        "5",
        "--wildcard-file",
        str(wc_linker1),
        str(ds_path),
    ]

    wc_linker2 = output_dir / "cutadapt_wc_L2.txt"
    log_linker2 = output_dir / "cutadapt_L2.log"
    dmuxL2 = [
        "cutadapt",
        "-g",
        f"linker2={linker2}",
        "-o",
        "/dev/null",
        "--action=lowercase",
        "--cores",
        str(cores),
        "--no-indels",
        "-e",
        "5",
        "--wildcard-file",
        str(wc_linker2),
        str(ds_path),
    ]

    logger.info("Starting dmuxL1")
    logger.debug("cutadapt L1 cmd: %s > %s", dmuxL1, log_linker1)
    with open(log_linker1, "w") as log_f:
        subprocess.run(dmuxL1, stdout=log_f, stderr=log_f, check=True)
    logger.info("Completed dmuxL1")

    logger.info("Starting dmuxL2")
    logger.debug("cutadapt L2 cmd: %s > %s", dmuxL2, log_linker2)
    with open(log_linker2, "w") as log_f:
        subprocess.run(dmuxL2, stdout=log_f, stderr=log_f, check=True)
    logger.info("Completed dmuxL2")

    return wc_linker1, log_linker1, wc_linker2, log_linker2


def write_onoff_table(onoff_df: pd.DataFrame, output_dir: Path) -> Path:
    onoff_path = output_dir / "onoff_tissue_table.csv"
    onoff_df.to_csv(onoff_path, index=False)
    return onoff_path
