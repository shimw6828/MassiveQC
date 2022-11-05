import os
from xopen import xopen
import logging
from pathlib import Path
from typing import Optional
from .fastq import Fastq, MixedUpReadsException, UnequalNumberReadsException
import pandas as pd
from .parser import remove_file

logger = logging.getLogger("MassiveQC")


class DownloadException(Exception):
    """Basic exception for problems downloading from SRA"""


class AbiException(Exception):
    """Basic exception when ABI file was downloaded from SRA"""


def check_and_compress_fastq(r1: str, QC_dir: str, r2: Optional[str] = None):
    """Check the reads quality and compress them in a new fastq file

    **parameter**

    r1: str
        The first fastq file

    QC_dir: str
        The directory of result fastq file after check_fq.

    r2: str or None
        The second fastq file

    **return**

    r1_gz: str

    r2_gz: str

    fq: Fastq class
    """
    fq = Fastq(r1, r2)
    r1_gz = os.path.join(QC_dir, os.path.basename(r1))
    if r2 is None:
        logger.info("Processing FASTQ as Single-End")
        run_as_se(fq, r1_gz)
        return r1_gz, fq
    else:
        logger.info("Processing FASTQ as Pair-End")
        r2_gz = os.path.join(QC_dir, os.path.basename(r2))
        run_as_pe(fq, r1_gz, r2_gz)
        return r1_gz, r2_gz, fq


def run_as_se(fq: Fastq, R1_out: str) -> None:
    """Check the reads quality and compress for single-end RNA-seq

    **parameter**

    fq: Fastq
        Fastq class

    R1_out: Path
        The result R1 fastq path

    **return**

    None

    """
    with xopen(R1_out, "wb") as file_out1:
        for read in fq.process():
            file_out1.write(read)
    if "abi_solid" in fq.flags:
        remove_file(R1_out)
        raise AbiException
    if "download_bad" in fq.flags:
        remove_file(R1_out)
        raise DownloadException("Empty FASTQ")
    if fq.libsize < 100_000:
        remove_file(R1_out)
        raise DownloadException("<100,000 reads")



def run_as_pe(fq: Fastq, R1_out: str, R2_out: str) -> None:
    """Check the reads quality and compress for single-end RNA-seq
    **parameter**

    fq: Fastq
        Fastq class

    R1_out: Path
        The result R1 fastq path

    R2_out: Path
        The result R2 fastq path

    **return**

    None
    """
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
        if "keep_R1" in fq.flags:
            run_as_se(fq, R1_out)
        elif "keep_R2" in fq.flags:
            run_as_se(fq, R2_out)


def save_output(feature_path, fq, SRR):
    """Save the output. This method will create a layout file contain 4 columns
    **parameter**

    feature_path: str
        The Feature dir path.

    fq: Fastq
        Fastq class.

    SRR: str
        SRR ID.

    **return**

    None
    """
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


def run_check_fq(SRR, SRA_path, QC_dir, feature_path):
    """Determine the layout of RNA-seq and Run check_and_compress_fastq
    **parameter**
    SRR: str
        SRR ID.
    SRA_path: str
        The path of fastq download directory.
    QC_dir: str
       The directory of result fastq file after check_fq.
    feature_path: str
        The directory of Feartures

    **return**
    None
    """
    file_list = [file for file in os.listdir(SRA_path) if file.startswith(SRR)]
    raw_fqs = []
    summary_file = os.path.join(feature_path, "layout", f"{SRR}.parquet")
    if len(file_list) == 2:
        logger.info("Pair-End QC")
        r1 = Path(SRA_path, f"{SRR}_1.fastq.gz")
        r2 = Path(SRA_path, f"{SRR}_2.fastq.gz")
        if os.path.exists(summary_file):
            return [r1.as_posix(), r2.as_posix()]
        r1_gz, r2_gz, fq = check_and_compress_fastq(r1=r1.as_posix(), QC_dir=QC_dir, r2=r2.as_posix())
        raw_fqs.append(r1)
        raw_fqs.append(r2)
    elif len(file_list) == 1:
        logger.info("Single-End QC")
        r1 = Path(SRA_path, f"{SRR}.fastq.gz")
        if os.path.exists(summary_file):
            return [r1]
        r1_gz, fq = check_and_compress_fastq(r1=r1.as_posix(), QC_dir=QC_dir)
        raw_fqs.append(r1)
    elif len(file_list) == 0:
        raise DownloadException
    save_output(feature_path, fq, SRR)
    return raw_fqs


def check_fq(SRR, SRA_path, QC_dir, feature_path):
    """This is the main function for checking the reads quality

    **parameter**

    SRR: str
        SRR ID.

    SRA_path: str
        The path of fastq download directory.

    QC_dir: str
       The directory of result fastq file after check_fq.

    feature_path: str
        The directory of Features.
    """
    try:
        raw_fqs = run_check_fq(SRR, SRA_path, QC_dir, feature_path)
        logger.info(f"Complete check {SRR} fastq file")
        return raw_fqs
    except AbiException:
        logger.warning(f"Flagging {SRR} as ABI Solid")
        raise AbiException
    except DownloadException as error:
        logger.warning(f"Flagging {SRR} as Download Bad: {error}")
        raise DownloadException
