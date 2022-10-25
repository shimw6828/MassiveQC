import argparse
import configparser
import logging
import os
import sys

import pandas as pd

logging.basicConfig(format='%(asctime)s %(message)s')
logger = logging.getLogger("MassiveQC")
logger.setLevel(logging.WARNING)
from pathlib import Path
from concurrent import futures
from tqdm import tqdm

sys.path.insert(0, "/home/mwshi/github/MassiveQC")
from MassiveQC.get_sra import get_sra
from MassiveQC.check_fq import check_fq
from MassiveQC.fastq_screen import fastq_screen
from MassiveQC.atropos import atropos
from MassiveQC.hisat2 import Hisat2
from MassiveQC.collectrnaseqmetrics import CollectRnaseqMetrics
from MassiveQC.markduplicates import MarkDuplicates
from MassiveQC.FeatureCounts import FeatureCounts
from MassiveQC.feature_store import check_done_sample, feature_store
from MassiveQC.detection import detection


def init_wd():
    Path(outdir).mkdir(exist_ok=True)
    Path(download_path).mkdir(exist_ok=True)
    feature_path.mkdir(exist_ok=True)
    (feature_path / "layout").mkdir(exist_ok=True)
    Path(QC_dir).mkdir(exist_ok=True)
    (feature_path / "fastq_screen").mkdir(exist_ok=True)
    (feature_path / "atropos").mkdir(exist_ok=True)
    Path(Bam_dir).mkdir(exist_ok=True)
    (feature_path / "hisat2").mkdir(exist_ok=True)
    (feature_path / "aln_stats").mkdir(exist_ok=True)
    (feature_path / "strand").mkdir(exist_ok=True)
    (feature_path / "rnaseqmetrics").mkdir(exist_ok=True)
    (feature_path / "genebody_coverage").mkdir(exist_ok=True)
    (feature_path / "markduplicates").mkdir(exist_ok=True)
    (feature_path / "count_summary").mkdir(exist_ok=True)
    Path(Count_dir).mkdir(exist_ok=True)
    (feature_path / "DoneSample").mkdir(exist_ok=True)


def process(SRR):
    logger.info(f"Start download {SRR}")
    fq_mode = get_sra(SRR, download_path, ascp_key)
    if only_download:
        return
    logger.info(f"Complete download {SRR}")
    logger.info(f"Start check {SRR} fastq file")
    # check_fq
    # Check if the result file exists.
    summary_file = feature_path / "layout" / f"{SRR}.parquet"
    if summary_file.exists():
        logger.info(f"{SRR} fastq file has been checked")
    else:
        check_fq(SRR, download_path, QC_dir, feature_path)
    logger.info(f"Complete check {SRR} fastq file")

    # fastq_screen
    fastq_screen_output = feature_path / "fastq_screen" / f"{SRR}.parquet"
    if fastq_screen_output.exists():
        logger.info(f"{SRR} fastq_screen step has been done")
    else:
        fastq_screen(SRR, QC_dir, feature_path.as_posix(), fastq_screen_config, THREADS)
    logger.info(f"Complete fastq_screen {SRR} fastq file")

    # atropos
    atropos_output = feature_path / "atropos" / f"{SRR}.parquet"
    if atropos_output.exists():
        logger.info(f"{SRR} atropos step has been done")
    else:
        atropos(feature_path.as_posix(), SRR, QC_dir, THREADS)

    # hisat2
    _hisat2 = Path(feature_path) / "hisat2" / f"{SRR}.parquet"
    _alnStat = Path(feature_path) / "aln_stats" / f"{SRR}.parquet"
    if _hisat2.exists() and _alnStat.exists():
        logger.info(f"{SRR} hisat2 step has been done")
    else:
        hisat_runner = Hisat2(feature_path.as_posix(), SRR, QC_dir, Bam_dir, THREADS, reference, splice=splice)
        hisat_runner.hisat2()

    # metrics
    strand = feature_path / "strand" / f"{SRR}.parquet"
    table = feature_path / "rnaseqmetrics" / f"{SRR}.parquet"
    coverage = feature_path / "genebody_coverage" / f"{SRR}.parquet"
    if strand.exists() and table.exists() and coverage.exists():
        logger.info(f"{SRR} collectrnaseqmetrics step has been done")
    else:
        metrics_runner = CollectRnaseqMetrics(feature_path.as_posix(), SRR, Bam_dir, THREADS, ref_flat, picard)
        metrics_runner.collectrnaseqmetrics()

    # markduplicates
    markdup = feature_path / "markduplicates" / f"{SRR}.parquet"
    if markdup.exists():
        logger.info(f"{SRR} markduplicates step has been done")
    else:
        markdup_runner = MarkDuplicates(feature_path.as_posix(), SRR, Bam_dir, THREADS, picard)
        markdup_runner.markduplicates()

    # FeatureCounts
    count_summary = feature_path / "count_summary" / f"{SRR}.parquet"
    if count_summary.exists():
        logger.info(f"{SRR} FeatureCounts step has been done")
    else:
        count_runner = FeatureCounts(feature_path.as_posix(), SRR, Bam_dir, Count_dir, gtf, THREADS)
        count_runner.FeatureCounts()

    # complete one srr, touch one file
    (feature_path / "DoneSample" / SRR).touch()
    return SRR


def local_thread(SRRs):
    init_wd()
    down_samples = os.listdir(feature_path / "DoneSample")
    pre_SRRs = [x for x in SRRs if x not in down_samples]
    with futures.ThreadPoolExecutor(max_workers=workers) as executor:
        tasks = [executor.submit(process, srr) for srr in pre_SRRs]
        for r in tqdm(futures.as_completed(tasks), total=len(tasks)):
            if r.exception():
                logger.error("One sample failed")
            else:
                res = r.result()



def get_arguments():
    # parse the config file
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument('-c', '--conf', type=str, required=False)
    pre_args, _ = pre.parse_known_args()
    config_args = {}
    if pre_args.conf:
        config = configparser.ConfigParser()
        with open(pre_args.conf, "r") as stream:
            config.read_string("[Defaults]\n" + stream.read())
        config_args.update(dict(config.items("Defaults")))
    parser = argparse.ArgumentParser(description='...', parents=[pre])
    parser.add_argument('-i', '--input', required=True, type=str, help='Input file, containing two columns srx and srr')
    parser.add_argument('-a', '--ascp_key', type=str, help='Locate aspera key. Default $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh')
    parser.add_argument('-f', '--fastq_screen_config', required=True, type=str, help="Path to the fastq_screen conf file, can be download from fastq_screen website")
    parser.add_argument('-g', '--gtf', required=True, type=str, help="Path to the GTF file with annotations")
    parser.add_argument('-x', '--ht2-idx', dest="ht2_idx", required=True, type=str, help="Hisat2 index filename prefix")
    parser.add_argument('-k', '--known-splicesite-infile', dest="known_splicesite_infile", type=str, help="Hisat2 splicesite file, provide a list of known splice sites")
    parser.add_argument('-p', '--picard', required=True, type=str, help="Path to picard.jar")
    parser.add_argument('-r', '--ref_flat', required=True, type=str, help="Path to refflat file")
    parser.add_argument('-o', '--outdir', required=True, type=str, help="Path to result output directory. If it doesn't exist, it will be created automatically")
    parser.add_argument('-w', '--workers', type=int, help="The number of simultaneous tasks", default=2)
    parser.add_argument('-t', '--THREADS', type=int, help="The number of threads for tools like Hisat2 in one task", default=4)
    parser.add_argument('-d', '--download', type=str, help="Path to SRA fastq files. The default is $outdir/download")
    parser.add_argument('--only_download', action="store_true", help="Only run the download step", default=False)
    parser.set_defaults(**config_args)
    if pre_args.conf:
        for action in parser._actions:
            if action.dest in config_args:
                action.required = False

    return parser.parse_args()


def main():
    args = get_arguments()
    global only_download
    only_download = args.only_download
    input_file = args.input.strip('"')
    global ascp_key
    ascp_key = args.ascp_key.strip('"')
    global gtf
    gtf = args.gtf.strip('"')
    global fastq_screen_config
    fastq_screen_config = args.fastq_screen_config.strip('"')
    global reference
    reference = args.ht2_idx.strip('"')
    global splice
    splice = args.known_splicesite_infile.strip('"')
    global picard
    picard = args.picard.strip('"')
    global ref_flat
    ref_flat = args.ref_flat.strip('"')
    global outdir
    outdir = args.outdir.strip('"')
    global workers
    workers = args.workers
    global THREADS
    THREADS = args.THREADS
    global download_path
    download_path = os.path.join(outdir, "download")
    global QC_dir
    QC_dir = os.path.join(outdir, "QC_dir")
    global Bam_dir
    Bam_dir = os.path.join(outdir, "Bam")
    global Count_dir
    Count_dir = os.path.join(outdir, "Count")
    global feature_path
    feature_path = os.path.join(outdir, "Features")
    feature_path = Path(feature_path)
    global done_sample
    done_sample = feature_path / "done_sample.txt"
    # init workshop
    init_wd()
    srr_df = pd.read_table(input_file)
    if len(srr_df.columns) == 1:
        # only have srr column
        srr_df.columns = ["srr"]
    elif len(srr_df.columns) == 2:
        # have srx and srr column
        srr_df.columns = ["srx", "srr"]

    SRRs = srr_df["srr"].values.tolist()
    logger.info(f"Start processing, {len(SRRs)} srrs will be processed")
    # run process local
    local_thread(SRRs)
    if not only_download:
        done_samples = check_done_sample(outdir)
        feature_store(done_samples, outdir)
        detection((feature_path / "features.parquet").as_posix())


if __name__ == "__main__":
    main()
