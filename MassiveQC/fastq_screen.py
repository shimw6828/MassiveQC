from pathlib import Path
import pandas as pd
import logging
from .command import run_command
from .parser import parse_fastq_screen

logger = logging.getLogger("MassiveQC")


def run_fastq_screen(config_file, feature_path, QC_dir, summary_file, SRR, THREADS=1):
    feature_screen = Path(feature_path) / "fastq_screen"
    layout_ = pd.read_parquet(summary_file).layout[0]
    QC_dir = Path(QC_dir)
    if layout_ == "PE" or layout_ == "keep_R1":
        fastq = QC_dir / f"{SRR}_1.fastq.gz"
        feature_file = feature_screen / f"{SRR}_1_screen.txt"
    elif layout_ == "keep_R2":
        fastq = QC_dir / f"{SRR}_2.fastq.gz"
        feature_file = feature_screen / f"{SRR}_1_screen.txt"
    else:
        fastq = QC_dir / f"{SRR}.fastq.gz"
        feature_file = feature_screen / f"{SRR}_screen.txt"

    screen(config_file, feature_screen, fastq, THREADS)
    output_file = feature_screen / f"{SRR}.parquet"
    summarize(feature_file, output_file, SRR)
    feature_file.unlink()
    (feature_screen / f"{feature_file.stem}.html").unlink()


def screen(config_file, feature_screen, fastq, THREADS: int) -> None:
    cmd = f"fastq_screen --outdir {feature_screen} " \
          f"--force --aligner bowtie2 --threads {THREADS} " \
          f"--conf {config_file} " \
          f"--subset 100000 " \
          f"{fastq} "
    logger.info(f"running {cmd}")
    log_info = run_command(cmd)
    if "Processing complete" not in log_info:
        logger.error(log_info)
        raise FastqScreenException


def summarize(feature_file: Path, output_file: Path, SRR: str) -> None:
    """Summarizes fastq screen results
        Calculates the number of reads mapping to each specific reference. Ignores reads the map to multiple references.
        Goes from this:
        | reference   |   multiple_hits_multiple_libraries_count |   multiple_hits_multiple_libraries_percent |   multiple_hits_one_library_count |   multiple_hits_one_library_percent |   one_hit_multiple_libraries_count |   one_hit_multiple_libraries_percent |   one_hit_one_library_count |   one_hit_one_library_percent |   reads_processed_count |   unmapped_count |   unmapped_percent |
        |:------------|-----------------------------------------:|-------------------------------------------:|----------------------------------:|------------------------------------:|-----------------------------------:|-------------------------------------:|----------------------------:|------------------------------:|------------------------:|-----------------:|-------------------:|
        | adapters    |                                       48 |                                       0.05 |                                 0 |                                0    |                                  0 |                                 0    |                           0 |                          0    |                   99973 |            99925 |              99.95 |
        | dm6         |                                     1713 |                                       1.71 |                              6278 |                                6.28 |                                224 |                                 0.22 |                       88393 |                         88.42 |                   99973 |             3365 |               3.37 |
        | ecoli       |                                        1 |                                       0    |                                 0 |                                0    |                                  0 |                                 0    |                           2 |                          0    |                   99973 |            99970 |             100    |
        ...
        To this:
        |            |   adapters_pct_reads_mapped |   dm6_pct_reads_mapped |   ecoli_pct_reads_mapped |   ercc_pct_reads_mapped |   hg19_pct_reads_mapped |   phix_pct_reads_mapped |   rRNA_pct_reads_mapped |   wolbachia_pct_reads_mapped |   yeast_pct_reads_mapped |
        |:-----------|----------------------------:|-----------------------:|-------------------------:|------------------------:|------------------------:|------------------------:|------------------------:|-----------------------------:|-------------------------:|
        | SRR0000001 |                           0 |                94.6966 |               0.00200054 |                       0 |               0.0160043 |                       0 |              0.00100027 |                            0 |               0.00500135 |
    """
    logger.info(F"Start extract {SRR} fastq screen result")
    df = parse_fastq_screen(feature_file).set_index("reference").fillna(0)
    summarized = (
        (
                (df.one_hit_one_library_count + df.multiple_hits_one_library_count)
                / df.reads_processed_count
                * 100
        )
            .rename(SRR)
            .rename_axis("")
            .to_frame()
            .T.rename_axis("srr")
    )
    summarized.columns = [f"{col}_pct_reads_mapped" for col in summarized.columns]
    summarized.to_parquet(output_file)
    logger.info("Complete, remove fastq screen result files")



class FastqScreenException(Exception):
    """Fastq Screen Processing Exception"""


def fastq_screen(SRR, QC_dir, feature_path, config_file, THREADS):
    """Use fastqscreen to identify RNA seq quality"""
    logger.info(f"Start fastq screen {SRR}")
    layout_file = Path(feature_path) / "layout" / f"{SRR}.parquet"
    try:
        run_fastq_screen(config_file, feature_path, QC_dir, layout_file, SRR, THREADS)
        logger.info(f"Complete fast_screen {SRR} fastq file")
    except FastqScreenException:
        logger.warning(f"{SRR}: fastq screen did not complete")
        raise FastqScreenException
