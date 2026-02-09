import logging
import sys

from datetime import datetime
from pathlib import Path

logger = logging.getLogger(__name__)


def format_wildcard_metrics(
    wc_name: str,
    total_reads: int,
    adapter_reads: int,
    unique_count: int,
    num_to_ninety: int,
    expected_count: int,
    whitelist_count: int,
    total_read_from_expected: float,
) -> str:
    return (
        f"\n######### Info for {wc_name}\n"
        f"Total Reads: {total_reads}  Reads with Adapter: {adapter_reads}\n"
        f"The number of unique strings in the 8mer column is {unique_count}.\n"
        f"Ninety percent (90%) of the reads come from total of {num_to_ninety} 8mers.\n"
        f"Total of {expected_count} out of {whitelist_count} expected 8-mers "
        f"accounted for {total_read_from_expected:.1%} of the reads\n"
    )


def format_hilo_metrics(
    label: str,
    total_hi_warn: int,
    total_lo_warn: int,
    total_mers: int,
) -> str:
    out = (
        f"\n######### Info for {label}\n"
        f"Total hiWarn: {total_hi_warn}\tTotal loWarn: {total_lo_warn}"
    )
    if total_mers > 0:
        out += (
            f"\nPct hiWarn: {(total_hi_warn / total_mers):.3f}\t"
            f"Pct loWarn: {(total_lo_warn / total_mers):.3f}"
        )
    else:
        out += "\nPct hiWarn: N/A\tPct loWarn: N/A"
    return out


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
