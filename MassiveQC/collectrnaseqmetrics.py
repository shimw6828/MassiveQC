import logging
import os.path
from pathlib import Path
import pandas as pd
from .command import run_command
from typing import Optional, Tuple
from .parser import parse_picardCollect_summary, parse_picardCollect_hist, remove_file

logger = logging.getLogger("MassiveQC")


class CollectRnaseqMetrics(object):
    """This class uses picard to extract RnaSeqMetrics information"""
    def __init__(self, feature_path: str, SRR: str, Bam_dir: str,
                 THREADS: int, ref_flat: str, picard: str, MEM: Optional[int] = 3):
        self.feature_path = Path(feature_path)
        self.layout = Path(feature_path) / f"{SRR}.parquet"
        self.SRR = SRR
        self.Bam_dir = Path(Bam_dir)
        self.THREADS = THREADS
        self.MEM = MEM
        self.ref_flat = ref_flat
        self.picard = picard

    def collectrnaseqmetrics(self):
        bam = self.Bam_dir / f"{self.SRR}.sorted.bam"
        unstranded, first, second = self.run_picard(bam)
        self.summarize(unstranded, first, second)
        remove_file(unstranded.as_posix())
        remove_file(first.as_posix())
        remove_file(second.as_posix())
        logger.info(f"{self.SRR} Complete metrics")

    def run_picard(self, bam: Path):
        cmd = (
            f"java -Xmx{self.MEM}g "
            f"-jar {self.picard} CollectRnaSeqMetrics "
            f"REF_FLAT={self.ref_flat} "
            f"INPUT={bam} OUTPUT={{out_file}} "
            f"STRAND={{strand}}"
        )
        logger.info(f"{self.SRR} Start CollectRnaSeqMetrics")
        # Unstranded
        un_out = self.feature_path / f"{self.SRR}.unstranded.txt"
        result_un = run_command(cmd.format(out_file=un_out, strand="NONE"))
        self._check_log(result_un, self.SRR)

        # First stranded
        first_out = self.feature_path / f"{self.SRR}.first_stranded.txt"
        result_first = run_command(cmd.format(out_file=first_out, strand="FIRST_READ_TRANSCRIPTION_STRAND"))
        self._check_log(result_first, self.SRR)

        # Second stranded
        second_out = self.feature_path / f"{self.SRR}.second_stranded.txt"
        result_second = run_command(cmd.format(out_file=second_out, strand="SECOND_READ_TRANSCRIPTION_STRAND"))
        self._check_log(result_second, self.SRR)
        return un_out, first_out, second_out

    def summarize(self, unstranded: Path, first: Path, second: Path):
        # Parse the flags
        if self._parse_stranded(first):
            strand_ = "same_strand"
        elif self._parse_stranded(second):
            strand_ = "opposite_strand"
        else:
            strand_ = "unstranded"

        idx = pd.Index([self.SRR], name="srr")
        # Write strand flag
        strand = self.feature_path / "strand" / f"{self.SRR}.parquet"
        strand_df = pd.DataFrame([[strand_]], index=idx, columns=["strand"])
        strand_df.to_parquet(strand)

        # Parse main table
        table = self.feature_path / "rnaseqmetrics" / f"{self.SRR}.parquet"
        table_df = self._parse_table(unstranded)
        table_df.index = idx
        table_df.to_parquet(table)

        # Parse genome coverage histogram
        coverage = self.feature_path / "genebody_coverage" / f"{self.SRR}.parquet"
        coverage_df = parse_picardCollect_hist(unstranded)
        coverage_df.index = idx
        coverage_df.to_parquet(coverage)

    @staticmethod
    def _parse_stranded(file_name: Path) -> bool:
        return (parse_picardCollect_summary(file_name).PCT_CORRECT_STRAND_READS >= 0.75)[0]

    @staticmethod
    def _parse_table(file_name: Path) -> pd.DataFrame:
        df = parse_picardCollect_summary(file_name)[
            [
                "PCT_CODING_BASES",
                "PCT_UTR_BASES",
                "PCT_INTRONIC_BASES",
                "PCT_INTERGENIC_BASES",
                "PCT_MRNA_BASES",
                "MEDIAN_CV_COVERAGE",
                "MEDIAN_5PRIME_BIAS",
                "MEDIAN_3PRIME_BIAS",
                "MEDIAN_5PRIME_TO_3PRIME_BIAS",
            ]
        ].fillna(0.0)
        df.PCT_CODING_BASES = df.PCT_CODING_BASES * 100
        df.PCT_UTR_BASES = df.PCT_UTR_BASES * 100
        df.PCT_INTRONIC_BASES = df.PCT_INTRONIC_BASES * 100
        df.PCT_INTERGENIC_BASES = df.PCT_INTERGENIC_BASES * 100
        df.PCT_MRNA_BASES = df.PCT_MRNA_BASES * 100

        df.columns = [x.lower() for x in df.columns]
        df.columns = [x.replace("pct_", "percent_") for x in df.columns]
        return df

    @staticmethod
    def _check_log(log: str, SRR) -> None:
        if "CollectRnaSeqMetrics done" not in log:
            raise PicardException(f"{SRR} not complete")


class PicardException(Exception):
    """Picard Processing Exception"""
