from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional

import barcodeqc.paths as paths
from barcodeqc.utils import require_executable, ExternalDependencyError

logger = logging.getLogger(__name__)

OUTPUT_DIR_SUFFIX = "_outputs"


def output_dir_from_sample_name(sample_name: str) -> Path:
    if sample_name.endswith(OUTPUT_DIR_SUFFIX):
        dir_name = sample_name
    else:
        dir_name = f"{sample_name}{OUTPUT_DIR_SUFFIX}"
    return Path.cwd() / dir_name


@dataclass(frozen=True)
class QCConfig:
    sample_name: str
    r2_path: Path
    barcode_set: Literal[
        "bc50", "bc96", "fg96", "bc220", "bc220_05-OCT", "bc220_20-MAY"
    ]
    sample_reads: int
    random_seed: int
    tissue_position_file: Optional[Path]
    output_dir: Path

    @classmethod
    def from_args(
        cls,
        sample_name: str,
        r2_path: Path,
        barcode_set: Literal[
            "bc50", "bc96", "fg96", "bc220", "bc220_05-OCT", "bc220_20-MAY"
        ],
        sample_reads: int,
        random_seed: int,
        tissue_position_file: Optional[Path],
    ) -> "QCConfig":
        output_dir = output_dir_from_sample_name(sample_name)
        if tissue_position_file is None:
            tissue_position_file = paths.BARCODE_PATHS[barcode_set]["positions"]
        return cls(
            sample_name=sample_name,
            r2_path=r2_path,
            barcode_set=barcode_set,
            sample_reads=sample_reads,
            random_seed=random_seed,
            tissue_position_file=tissue_position_file,
            output_dir=output_dir,
        )

    def validate(self) -> None:
        if not self.r2_path.exists():
            raise FileNotFoundError(
                f"fastq file path does not exist: {self.r2_path}"
            )
        if self.tissue_position_file is not None:
            if not self.tissue_position_file.exists():
                raise FileNotFoundError(
                    f"Could not find tissue_postion file: {self.tissue_position_file}"
                )

        # Ensure seqtk installed in PATH
        seqtk = require_executable("seqtk")
        logger.debug(f"Using {seqtk} for subsampling.")
