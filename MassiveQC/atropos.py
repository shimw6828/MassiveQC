from pathlib import Path
import pandas as pd
import logging
from typing import Tuple
import re
from .command import run_command

logger = logging.getLogger("MassiveQC")


def atropos(feature_path: str, SRR: str, QC_dir: str, THREADS: int):
    """This function is main part which is used to trim the reads.
    Filter reads that are less than 25bp.

    **parameter**
    feature_path: str
        The directory of Features.
    SRR: str
        SRR ID.
    QC_dir: str
       The directory of result fastq file after check_fq.
    THREADS: int
        Thread number for atropos.
    **return**
    None
    """
    logger.info(f"Start atropos {SRR}")
    try:
        QC_dir = Path(QC_dir)
        layout = Path(feature_path) / "layout" / f"{SRR}.parquet"
        layout_ = pd.read_parquet(layout).layout[0]
        results = run_atropos(layout_, SRR, QC_dir, THREADS)
        output = Path(feature_path) / "atropos" / f"{SRR}.parquet"
        summarize(results, output, SRR)
        if layout_ == "PE":
            r1_trim = QC_dir / f"{SRR}_1.trim.fastq.gz"
            r2_trim = QC_dir / f"{SRR}_2.trim.fastq.gz"
            if r1_trim.exists() and r2_trim.exists():
                (QC_dir / f"{SRR}_1.fastq.gz").unlink()
                (QC_dir / f"{SRR}_2.fastq.gz").unlink()
        elif layout_ == "Keep_R1":
            r1_trim = QC_dir / f"{SRR}_1.trim.fastq.gz"
            if r1_trim.exists():
                (QC_dir / f"{SRR}_1.fastq.gz").unlink()
        elif layout_ == "Keep_R2":
            r1_trim = QC_dir / f"{SRR}_2.trim.fastq.gz"
            if r1_trim.exists():
                (QC_dir / f"{SRR}_2.fastq.gz").unlink()
        else:
            r1_trim = QC_dir / f"{SRR}.trim.fastq.gz"
            if r1_trim.exists():
                (QC_dir / f"{SRR}.fastq.gz").unlink()
        logger.info(f"Complete atropos {SRR}")
    except AtroposException as error:
        logger.warning(f"Flagging {SRR} as Atropos Bad")
        raise AtroposException("Atropos Bad")


def run_atropos(layout_, SRR, QC_dir: Path, THREADS) -> str:
    """Run the atropos command in shell"""
    # adapters = os.path.join(os.path.abspath(__file__), "sequencing_adapters.fa")
    if layout_ == "PE":
        r1 = QC_dir / f"{SRR}_1.fastq.gz"
        r2 = QC_dir / f"{SRR}_2.fastq.gz"
        r1_trim = QC_dir / f"{SRR}_1.trim.fastq.gz"
        r2_trim = QC_dir / f"{SRR}_2.trim.fastq.gz"
        cmd = f"atropos trim " \
              f"-q 20 --minimum-length 25 " \
              f"--threads {THREADS} " \
              f"-pe1 {r1} -pe2 {r2} -o {r1_trim} -p {r2_trim}"
    elif layout_ == "Keep_R1":
        r1 = QC_dir / f"{SRR}_1.fastq.gz"
        r1_trim = QC_dir / f"{SRR}_1.trim.fastq.gz"
        cmd = f"atropos trim " \
              "-q 20 -U 0 --minimum-length 25 " \
              f"--threads {THREADS} " \
              f"-se {r1} -o {r1_trim}"
    elif layout_ == "Keep_R2":
        r2 = QC_dir / f"{SRR}_2.fastq.gz"
        r2_trim = QC_dir / f"{SRR}_2.trim.fastq.gz"
        cmd = f"atropos trim " \
              "-q 20 -U 0 --minimum-length 25 " \
              f"--threads {THREADS} " \
              f"-se {r2} -o {r2_trim}"
    else:
        r1 = QC_dir / f"{SRR}.fastq.gz"
        r1_trim = QC_dir / f"{SRR}.trim.fastq.gz"
        cmd = f"atropos trim " \
              "-q 20 -U 0 --minimum-length 25 " \
              f"--threads {THREADS} " \
              f"-se {r1} -o {r1_trim}"
    logger.info(f"running {cmd}")
    results = run_command(cmd, verbose=True)
    return results


def summarize(log: str, output: Path, SRR: str) -> None:
    """Extract features from atropos result"""
    df = pd.DataFrame(
        [parse_atropos_log(log, SRR)],
        index=pd.Index([SRR], name="srr"),
        columns=["total_processed", "total_written", "too_short"],
    )
    if df.total_written[0] < 1_000:
        raise AtroposException("<1,000 reads")
    df.to_parquet(output)


def parse_atropos_log(log_text: str, SRR: str) -> Tuple[int, int, int]:
    """Parse atropos log information
       Example Atropos output
        # Example of SE
        Reads                                 records  fraction
        ----------------------------------- --------- ---------
        Total reads processed:                 85,319
        Reads with adapter:                     4,725      5.5%
        Reads that were too short:              1,879      2.2%
        Reads written (passing filters):       83,440     97.8%
        # Example of PE
        Pairs                                records fraction
        ----------------------------------- -------- --------
        Total read pairs processed:            1,354
        Read 1 with adapter:                   106     7.8%
        Read 2 with adapter:                   123     9.1%
        Pairs that were too short:                16     1.2%
        Pairs written (passing filters):       1,338    98.8%
    """
    if "ERROR" in log_text:
        logger.error(log_text)
        logger.warning(f"{SRR}: Atropos reported an error")
    try:
        log_text = log_text.replace(",", "")
        tot_processed = int(re.findall(r"Total read.*processed:\s+(\d+)", log_text)[0])
        tot_written = int(
            re.findall(r"[Read|Pair]s written \(passing filters\):\s+(\d+)", log_text)[0]
        )
        too_short = int(re.findall(r"[Read|Pair]s that were too short:\s+(\d+)", log_text)[0])
        return tot_processed, tot_written, too_short
    except IndexError:
        raise AtroposException("Unable to parse log")


class AtroposException(Exception):
    """Basic Atropos Exception"""
