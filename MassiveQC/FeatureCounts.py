import logging
from pathlib import Path
import pandas as pd
from .command import run_command
from .parser import remove_file
import numpy as np

logger = logging.getLogger("MassiveQC")


class FeatureCounts(object):
    """Run FeatureCounts and extract count features"""
    def __init__(self, feature_path: str, SRR: str, Bam_dir: str, Count_dir: str,
                 gtf: str, THREADS: int):
        self.feature_path = Path(feature_path)
        self.layout = self.feature_path / "layout" / f"{SRR}.parquet"
        self.strand = self.feature_path / "strand" / f"{SRR}.parquet"
        self.SRR = SRR
        self.Bam_dir = Path(Bam_dir)
        self.Count_dir = Path(Count_dir)
        self.gtf = gtf
        self.THREADS = THREADS

    def FeatureCounts(self):
        counts, jcounts = self.run_featureCounts()
        self.summarize(counts, jcounts)
        return (self.Bam_dir / f"{self.SRR}.sorted.bam").as_posix()

    def run_featureCounts(self):
        # Look up Layout
        layout_ = pd.read_parquet(self.layout).layout[0]
        if layout_ == "PE":
            params = "-p -P -C -J -B "
        else:
            params = "-J "

        # Look up strand
        strand_ = pd.read_parquet(self.strand).strand[0]
        if strand_ == "same_strand":
            params += "-s 1"
        elif strand_ == "opposite_strand":
            params += "-s 2"
        else:
            params += "-s 0"

        counts = self.Count_dir / f"{self.SRR}.counts"
        jcounts = self.Count_dir / f"{self.SRR}.counts.jcounts"
        summary = self.Count_dir / f"{self.SRR}.counts.summary"
        bam = self.Bam_dir / f"{self.SRR}.sorted.bam"

        cmd = (
            f"featureCounts  -T {self.THREADS} {params} "
            f"-a {self.gtf} -o {counts} "
            f"{bam}"
        )
        logger.info(f"{cmd}")
        logger.info(f"{self.SRR} start featureCounts")
        _result = run_command(cmd)
        self._check_log(_result, self.SRR)
        remove_file(summary.as_posix())
        logger.info(f"{self.SRR} complete featureCounts")
        return counts, jcounts

    def summarize(self, counts: Path, jcounts: Path):
        gene_counts = self._get_counts(counts)
        genic_reads = gene_counts.sum()
        percent_genes_on = (gene_counts > 0).mean() * 100

        junction_counts = self._get_counts(jcounts)
        junction_reads = junction_counts.sum()
        number_junctions_on = junction_counts.shape[0]

        df = pd.DataFrame(
            [[genic_reads, percent_genes_on, junction_reads, number_junctions_on]],
            columns=[
                "number_genic_reads",
                "percent_genes_on",
                "number_junction_reads",
                "number_junctions_on",
            ],
            index=pd.Index([self.SRR], name="srr"),
        )
        output_file = self.feature_path / "count_summary" / f"{self.SRR}.parquet"
        df.to_parquet(output_file)

    @staticmethod
    def _get_counts(file_name: Path) -> pd.DataFrame:
        col_name = pd.read_table(file_name, comment="#", nrows=1).columns[-1]
        return pd.read_table(file_name, comment="#", usecols=[col_name], dtype=np.int64).values

    @staticmethod
    def _check_log(log: str, SRR) -> None:
        if not (
                "Read assignment finished" in log or "Summary of counting results" in log
        ):
            raise FeatureCountsException(f"{SRR} not complete")


class FeatureCountsException(Exception):
    """Fastq Screen Processing Exception"""
