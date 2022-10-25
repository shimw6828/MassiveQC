import os, argparse
import sys
from xopen import xopen
import logging
from pathlib import Path
from typing import Optional
from .fastq import Fastq, MixedUpReadsException, UnequalNumberReadsException
import pandas as pd

logger = logging.getLogger("MassiveQC")


class DownloadException(Exception):
    """Basic exception for problems downloading from SRA"""


class AbiException(Exception):
    """Basic exception when ABI file was downloaded from SRA"""


def check_and_compress_fastq(r1: str, QC_dir: str, r2: Optional[str] = None):
    fq = Fastq(r1, r2)
    r1_gz = os.path.join(QC_dir, os.path.basename(r1))
    if r2 is None:
        # r2_gz = QC_dir / "empty_"+os.path.basename(r2)
        logger.info("Processing FASTQ as Single-End")
        run_as_se(fq, r1_gz)
        return r1_gz, fq
    else:
        logger.info("Processing FASTQ as Pair-End")
        r2_gz = os.path.join(QC_dir, os.path.basename(r2))
        run_as_pe(fq, r1_gz, r2_gz)
        return r1_gz, r2_gz, fq


def run_as_se(fq: Fastq, R1_out: Path) -> None:
    with xopen(R1_out, "wb") as file_out1:
        for read in fq.process():
            file_out1.write(read)
    if "abi_solid" in fq.flags:
        raise AbiException
    if "download_bad" in fq.flags:
        raise DownloadException("Empty FASTQ")
    if fq.libsize < 100_000:
        raise DownloadException("<100,000 reads")


def run_as_pe(fq: Fastq, R1_out: Path, R2_out: Path) -> None:
    try:
        with xopen(R1_out, "wb") as file_out1, xopen(R2_out, "wb") as file_out2:
            for read1, read2 in fq.process():
                file_out1.write(read1)
                file_out2.write(read2)
        if "abi_solid" in fq.flags:
            raise AbiException
        if "download_bad" in fq.flags:
            raise DownloadException("Empty FASTQ")
        if fq.libsize < 100_000:
            raise DownloadException("<100,000 reads")
    except (UnequalNumberReadsException, MixedUpReadsException):
        remove_file(R1_out)
        remove_file(R2_out)
        run_as_se(fq, R1_out, R2_out)


def save_output(feature_path, fq, SRR):
    summary_file = os.path.join(feature_path, "layout", f"{SRR}.parquet")
    layout = fq.flags.intersection(set(["SE", "PE", "keep_R1", "keep_R2"])).pop()
    idx = pd.Index([SRR], name="srr")
    if isinstance(fq.avgReadLen, list):
        r1, r2 = fq.avgReadLen
    else:
        r1, r2 = fq.avgReadLen, 0.0
    df = pd.DataFrame([[layout, fq.libsize, r1, r2]], index=[idx],
                      columns=["layout", "libsize", "avgLen_R1", "avgLen_R2"])
    df.to_parquet(summary_file)


def remove_file(file_name: str):
    if file_name is None:
        return
    pth = Path(file_name)
    if pth.exists() & pth.is_file():
        pth.unlink()


def remove_outputs(outputs) -> None:
    for output in outputs:
        remove_file(output)


def main(SRR, SRA_path, QC_dir, feature_path):
    file_list = [file for file in os.listdir(SRA_path) if file.startswith(SRR)]
    if len(file_list) == 2:
        logger.info("Pair-End QC")
        r1 = Path(SRA_path, f"{SRR}_1.fastq.gz")
        r2 = Path(SRA_path, f"{SRR}_2.fastq.gz")
        r1_gz, r2_gz, fq = check_and_compress_fastq(r1=r1.as_posix(), QC_dir=QC_dir, r2=r2.as_posix())
    elif len(file_list) == 1:
        logger.info("Single-End QC")
        r1 = Path(SRA_path, f"{SRR}.fastq.gz")
        r1_gz, fq = check_and_compress_fastq(r1=r1.as_posix(), QC_dir=QC_dir)
    elif len(file_list) == 0:
        raise DownloadException
    save_output(feature_path, fq, SRR)


def check_fq(SRR, SRA_path, QC_dir, feature_path):
    # bad_summary_file = os.path.join(feature_path, "layout", f"bad_{SRR}.parquet")
    try:
        main(SRR, SRA_path, QC_dir, feature_path)
    except AbiException:
        logger.warning(f"Flagging {SRR} as ABI Solid")
        # idx = pd.Index([SRR], name="srr")
        # df = pd.DataFrame([["Abi", 0, 0, 0]], index=[idx],
        #                   columns=["layout", "libsize", "avgLen_R1", "avgLen_R2"])
        # df.to_parquet(bad_summary_file)
        raise AbiException
    except DownloadException as error:
        logger.warning(f"Flagging {SRR} as Download Bad: {error}")
        # idx = pd.Index([SRR], name="srr")
        # df = pd.DataFrame([["Download_bad", 0, 0, 0]], index=[idx],
        #                   columns=["layout", "libsize", "avgLen_R1", "avgLen_R2"])
        # df.to_parquet(bad_summary_file)
        raise DownloadException
