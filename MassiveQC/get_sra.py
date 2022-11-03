"""Extract SRA file into Gziped FASTQs"""
import logging
import os
from pathlib import Path
from pysradb.sraweb import SRAweb
from .command import run_command

db = SRAweb()
logger = logging.getLogger("MassiveQC")


class DownloadException(Exception):
    """Basic exception for problems downloading from SRA"""


def _find_aspera_keypath(ascp_key=None):
    """Locate aspera key.
    Parameters
    ----------
    aspera_dir: str
                Location to aspera directory (optional)
    Returns
    -------
    aspera_keypath: str
                    Location to aspera key
    """
    if ascp_key is None:
        aspera_dir = os.path.join(os.path.expanduser("~"), ".aspera")
        aspera_keypath = os.path.join(
            aspera_dir, "connect", "etc", "asperaweb_id_dsa.openssh"
        )
    else:
        aspera_keypath = ascp_key
    if os.path.isfile(aspera_keypath):
        return aspera_keypath


def verify_sra_download(log, SRR):
    if "error" in log.lower():
        logger.warning("Warning: ascp download failed, start download with wget")
        raise DownloadException("Download Failed")
    if "failed" in log.lower() or "unable" in log.lower():
        logger.error("wget download also failed")
        raise DownloadException(f"{SRR} Download Failed")
    if "command not found" in log.lower():
        raise DownloadException("Failed ascp not install")


def sra_ascp(SRR: str, download_path: str, ascp_key) -> Path:
    ascp = "ascp -k1 -T -l 300m -P33001 -i"
    record = db.sra_metadata(SRR, detailed=True)
    record = record.dropna(axis=1, how="any")
    for a in record.to_dict('records'):
        if SRR == a["run_accession"]:
            record = a
    ena_cols = [x for x in list(record.keys()) if "ena_fastq_ftp" in x]
    fastq_col = [x for x in list(record.keys()) if "ena_fastq_http" in x]
    if len(ena_cols) == 0:
        logger.warning(f"{SRR} can not find fastq file on EBI")
        raise DownloadException(f"{SRR} can not find fastq file on EBI")
    filenum = 0
    for ena, fastq in zip(ena_cols, fastq_col):
        try:
            download_url = record[ena]
            if download_url == "":
                continue
            if filenum == 0:
                logger.info(f"{SRR} first file start download")
            else:
                logger.info(f"{SRR} second file start download")
            cmd = "{} {} {} {}".format(
                ascp, _find_aspera_keypath(ascp_key), download_url, download_path
            )
            logger.info(f"running {cmd}")
            log = run_command(cmd)
            verify_sra_download(log, SRR)
        except:
            download_url = record[fastq]
            if download_url == "":
                continue
            if filenum == 0:
                logger.info(f"{SRR} first file start download")
            else:
                logger.info(f"{SRR} second file start download")
            cmd = "wget -N -c -q --timeout=120 -P {} {}".format(
                download_path, download_url
            )
            logger.info(f"running {cmd}")
            log = run_command(cmd)
            verify_sra_download(log, SRR)
        filenum += 1


def get_sra(SRR: str, download_path: str, ascp_key=None) -> Path:
    """Download sra fastq

    **parameter**
    SRR: str
        SRR ID.
    download_path: str
        Download directory
    ascp_key: str
        Location to aspera directory (optional)

    **return**
    Path
    """
    pe_r1 = Path(download_path) / f"{SRR}_1.fastq.gz"
    pe_r2 = Path(download_path) / f"{SRR}_2.fastq.gz"
    se_r1 = Path(download_path) / f"{SRR}.fastq.gz"
    if Path(se_r1).exists():
        logger.info("The file already exists")
        return
    if Path(pe_r1).exists() and Path(pe_r2).exists():
        logger.info("The files already exists")
        return
    try:
        sra_ascp(SRR, download_path, ascp_key)
    except:
        srrs = os.listdir(download_path)
        for srr_file in srrs:
            if SRR in srr_file:
                os.unlink(os.path.join(download_path, srr_file))
        raise DownloadException(f"{SRR} download failed")
