import logging
from pathlib import Path
import pandas as pd
from .command import run_command
from typing import Optional, Tuple
from .parser import parse_hisat2, parse_samtools_stats, parse_bamtools_stats, remove_file

logger = logging.getLogger("MassiveQC")


class Hisat2(object):
    """This class can run hisat2 and extract the alignment summary"""
    def __init__(self, feature_path: str, SRR: str, QC_dir: str, Bam_dir: str,
                 THREADS: int, reference: str, strand: Optional[str] = None,
                 splice: Optional[str] = None):
        self.feature_path = Path(feature_path)
        self.layout = Path(feature_path) / "layout" / f"{SRR}.parquet"
        self.SRR = SRR
        self.QC_dir = Path(QC_dir)
        self.Bam_dir = Path(Bam_dir)
        self.THREADS = THREADS
        self.reference = reference
        self.strand = strand
        self.splice = splice
        self.r1 = None
        self.r2 = None

    def hisat2(self):
        logger.info(f"")
        layout_ = pd.read_parquet(self.layout).layout[0]
        trim_fqs = []
        if layout_ == "PE":
            self.r1 = self.QC_dir / f"{self.SRR}_1.trim.fastq.gz"
            self.r2 = self.QC_dir / f"{self.SRR}_2.trim.fastq.gz"
        elif layout_ == "Keep_R1":
            self.r1 = self.QC_dir / f"{self.SRR}_1.trim.fastq.gz"
            r2 = None
        elif layout_ == "Keep_R2":
            self.r1 = self.QC_dir / f"{self.SRR}_2.trim.fastq.gz"
            self.r2 = None
        else:
            self.r1 = self.QC_dir / f"{self.SRR}.trim.fastq.gz"
            self.r2 = None
        results, sam = self.run_hisat2()
        # return results
        _hisat2 = self.feature_path / "hisat2" / f"{self.SRR}.parquet"
        self.check_hisat(results, _hisat2, self.SRR)
        bam, bai = self.compress_sort_and_index(sam)
        remove_file(sam.as_posix())
        _alnStat = self.feature_path / "aln_stats" / f"{self.SRR}.parquet"
        self.alignment_stats(bam, _alnStat)
        trim_fqs.append(self.r1.as_posix())
        trim_fqs.append(self.r2.as_posix())
        return trim_fqs

    def run_hisat2(self):
        if self.r2:
            fastqs = f"-1 {self.r1} -2 {self.r2}"
        else:
            fastqs = f"-U {self.r1}"
        # Look up strand if it is there
        strand_param = ""
        if self.strand:
            strand_ = pd.read_parquet(self.strand).strand[0]
            if (self.layout_ == "PE") & (strand_ == "first_strand"):
                strand_param = "--rna-strandness FR"
            elif (self.layout_ == "PE") & (strand_ == "second_strand"):
                strand_param = "--rna-strandness RF"
            elif strand_ == "first_strand":
                strand_param = "--rna-strandness F"
            elif strand_ == "second_strand":
                strand_param = "--rna-strandness R"
        splice_param = ""
        if self.splice:
            splice_param = f"--known-splicesite-infile {self.splice}"
        sam = self.Bam_dir / f"{self.SRR}.sam"
        cmd = (
            "hisat2 "
            f"-x {self.reference} "
            f"{fastqs} "
            f"--threads {self.THREADS} "
            "--max-intronlen 300000 "
            f"{strand_param} "
            f"{splice_param} "
            f"-S {sam} "
        )
        logger.info(f"{self.SRR} Start Hisat2 alignment")
        logger.info(f"{cmd}")
        results = run_command(cmd)
        logger.info(f"{self.SRR} Complete Hisat2 alignment")
        return results, sam

    def compress_sort_and_index(self, sam: Path) -> Tuple[Path, Path]:
        sorted_bam = self.Bam_dir / f"{self.SRR}.sorted.bam"
        sorted_bai = self.Bam_dir / f"{self.SRR}.sorted.bam.bai"
        cmd = f"samtools view -Sb -q 20 --threads {self.THREADS} {sam} " \
              f"|samtools sort -l 9 --output-fmt BAM " \
              f" --threads {self.THREADS} -o {sorted_bam} && " \
              f"samtools index {sorted_bam}"
        logger.info(cmd)
        results = run_command(cmd)
        if "error" in results.lower():
            logger.warning(f"{self.SRR} samtoots error")
            raise Exception(f"samtoots error")
        else:
            return sorted_bam, sorted_bai

    def alignment_stats(self, bam: Path, output_file: Path):
        logger.info(f"samtools stats {bam}")
        result1 = run_command(f"samtools stats {bam}", verbose=False)
        logger.info(f"bamtools stats -in {bam}")
        result2 = run_command(f"bamtools stats -in {bam}", verbose=False)
        # Summarize
        df = pd.concat([self._samtools(result1), self._bamtools(result2)], axis=1, sort=False)
        df.index = pd.Index([self.SRR])
        df.to_parquet(output_file)

    @staticmethod
    def check_hisat(log: str, output_file: str, srr) -> None:
        """Example Hisat2 output
            # Example SE
            83440 reads; of these:
            83440 (100.00%) were unpaired; of these:
                3605 (4.32%) aligned 0 times
                76302 (91.45%) aligned exactly 1 time
                3533 (4.23%) aligned >1 times
            95.68% overall alignment rate
            # Example PE
            1338 reads; of these:
            1338 (100.00%) were paired; of these:
                468 (34.98%) aligned concordantly 0 times
                839 (62.71%) aligned concordantly exactly 1 time
                31 (2.32%) aligned concordantly >1 times
                ----
                468 pairs aligned concordantly 0 times; of these:
                22 (4.70%) aligned discordantly 1 time
                ----
                446 pairs aligned 0 times concordantly or discordantly; of these:
                892 mates make up the pairs; of these:
                    807 (90.47%) aligned 0 times
                    74 (8.30%) aligned exactly 1 time
                    11 (1.23%) aligned >1 times
            69.84% overall alignment rate
        """
        df = parse_hisat2(log).fillna(0)
        df.index = pd.Index([srr])
        df.to_parquet(output_file)
        uniquely_aligned = (
                df.iloc[0, :]["num_concordant_reads_uniquely_aligned"]
                + df.iloc[0, :]["num_uniquely_aligned"]
        )

        per_aligned = df.iloc[0, :]["per_alignment"]

        if (per_aligned < 1) | (uniquely_aligned < 1000):
            raise Exception(f"Poor alignment: {uniquely_aligned:,} ({per_aligned}%)")

    @staticmethod
    def _samtools(stats_file: Path) -> pd.DataFrame:
        return parse_samtools_stats(stats_file)[
            [
                "reads_MQ0",
                "average_quality",
                "insert_size_average",
                "insert_size_standard_deviation",
                "inward_oriented_pairs",
                "outward_oriented_pairs",
                "pairs_with_other_orientation",
                "pairs_on_different_chromosomes",
            ]
        ]

    @staticmethod
    def _bamtools(stats_file: Path) -> pd.DataFrame:
        return parse_bamtools_stats(stats_file)[["Percent Forward", "Percent Reverse"]]
