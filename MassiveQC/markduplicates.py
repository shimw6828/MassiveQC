import logging
from pathlib import Path
import pandas as pd
from .command import run_command
from typing import Optional
from .parser import parse_picard_markduplicate_metrics, remove_file
import numpy as np

logger = logging.getLogger("MassiveQC")


class MarkDuplicates(object):
    """run picard MarkDuplicates and extract the result"""
    def __init__(self, feature_path: str, SRR: str, Bam_dir: str,
                 THREADS: int, picard: str, MEM: Optional[int] = 3):
        self.feature_path = Path(feature_path)
        self.layout = Path(feature_path) / f"{SRR}.parquet"
        self.SRR = SRR
        self.Bam_dir = Path(Bam_dir)
        self.THREADS = THREADS
        self.MEM = MEM
        self.picard = picard

    def markduplicates(self):
        bam = self.Bam_dir / f"{self.SRR}.sorted.bam"
        metrics = self.run_markduplicates(bam)
        self.summarize(metrics)

    def run_markduplicates(self, bam: Path):
        dedup_bam = self.Bam_dir / f"{self.SRR}.dedup.bam"
        metrics = self.feature_path / f"{self.SRR}.metrics"
        cmd = (
            f"java -Xmx{self.MEM}g "
            f"-jar {self.picard} MarkDuplicates "
            f"INPUT={bam} "
            f"OUTPUT={dedup_bam} "
            f"METRICS_FILE={metrics}"
        )
        logger.info(f"{cmd}")
        logger.info(f"{self.SRR} Start markduplicates")
        _result = run_command(cmd)
        self._check_log(_result, self.SRR)
        remove_file(dedup_bam.as_posix())
        return metrics

    def summarize(self, metrics: Path):
        dtypes = {
            "UNPAIRED_READS_EXAMINED": np.int64,
            "READ_PAIRS_EXAMINED": np.int64,
            "UNPAIRED_READ_DUPLICATES": np.int64,
            "READ_PAIR_DUPLICATES": np.int64,
            "PERCENT_DUPLICATION": np.float64,
            "ESTIMATED_LIBRARY_SIZE": np.int64,
        }
        df = parse_picard_markduplicate_metrics(metrics)[dtypes.keys()].fillna(0).astype(dtypes)
        df.PERCENT_DUPLICATION = df.PERCENT_DUPLICATION * 100
        df.columns = [col.lower() for col in df.columns]
        df.index = pd.Index([self.SRR], name="srr")
        output_file = self.feature_path / "markduplicates" / f"{self.SRR}.parquet"
        remove_file(metrics.as_posix())
        df.to_parquet(output_file)

    @staticmethod
    def _check_log(log: str, SRR) -> None:
        if "MarkDuplicates done" not in log:
            raise PicardException(f"{SRR} not complete")


class PicardException(Exception):
    """Picard Processing Exception"""
