import argparse
import configparser
import logging
import os

logging.basicConfig(format='%(asctime)s %(message)s')
logger = logging.getLogger("MassiveQC")
logger.setLevel(logging.INFO)
from pathlib import Path

# sys.path.insert(0, "/home/mwshi/github/MassiveQC")
from .get_sra import get_sra
from .check_fq import check_fq
from .fastq_screen import fastq_screen
from .atropos import atropos
from .hisat2 import Hisat2
from .collectrnaseqmetrics import CollectRnaseqMetrics
from .markduplicates import MarkDuplicates
from .FeatureCounts import FeatureCounts
from .parser import remove_file


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
    if not skip_download:
        fq_mode = get_sra(SRR, download_path, ascp_key)
        logger.info(f"Complete download {SRR}")
    if only_download:
        return

    # check_fq
    # Check if the result file exists.
    logger.info(f"Start check {SRR} fastq file")
    summary_file = feature_path / "layout" / f"{SRR}.parquet"
    if summary_file.exists():
        logger.info(f"{SRR} fastq file has been checked")
    else:
        raw_fqs = check_fq(SRR, download_path, QC_dir, feature_path)
        # remove the raw fastq
        if remove_fastq:
            for k in raw_fqs:
                logger.info(f"Remove {k}")
                remove_file(k)

    # fastq_screen
    fastq_screen_output = feature_path / "fastq_screen" / f"{SRR}.parquet"
    if fastq_screen_output.exists():
        logger.info(f"{SRR} fastq_screen step has been done")
    else:
        fastq_screen(SRR, QC_dir, feature_path.as_posix(), fastq_screen_config, THREADS)

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
        trim_fqs = hisat_runner.hisat2()
        if remove_fastq:
            for k in trim_fqs:
                logger.info(f"Remove {k}")
                remove_file(k)

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
        bam_file = count_runner.FeatureCounts()
        if remove_bam:
            remove_file(bam_file)

    # complete one srr, touch one file
    (feature_path / "DoneSample" / SRR).touch()
    return SRR


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
    parser = argparse.ArgumentParser(description='...', formatter_class=argparse.RawTextHelpFormatter, parents=[pre])

    parser.add_argument('-s', '--srr', required=True, type=str, help='SRR id')
    parser.add_argument('-a', '--ascp_key', type=str, help='Locate aspera key. Default=$HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh')
    parser.add_argument('-f', '--fastq_screen_config', required=True, type=str, help="Path to the fastq_screen conf file, can be download from fastq_screen website")
    parser.add_argument('-g', '--gtf', required=True, type=str, help="Path to the GTF file with annotations")
    parser.add_argument('-x', '--ht2-idx', dest="ht2_idx", required=True, type=str, help="Hisat2 index filename prefix")
    parser.add_argument('-k', '--known-splicesite-infile', dest="known_splicesite_infile", type=str, help="Hisat2 splicesite file, provide a list of known splice sites")
    parser.add_argument('-p', '--picard', required=True, type=str, help="Path to picard.jar")
    parser.add_argument('-r', '--ref_flat', required=True, type=str, help="Path to refflat file")
    parser.add_argument('-o', '--outdir', required=True, type=str, help="Path to result output directory. If it doesn't exist, it will be created automatically")
    parser.add_argument('-t', '--THREADS', type=int, help="The number of threads for tools like Hisat2 in one task", default=4)
    parser.add_argument('-d', '--download', type=str, help="Path to SRA fastq files. The default is $OUTDIR/download")
    parser.add_argument('--only_download', action="store_true", help="Only run the download step", default=False)
    parser.add_argument('--skip_download', action="store_true", help="Skip the download step", default=False)
    parser.add_argument('--remove_fastq', action="store_true", help="Don't remain the fastq after running hisat2", default=False)
    parser.add_argument('--remove_bam', action="store_true", help="Don't remain the bam after running FeatureCounts", default=False)
    parser.set_defaults(**config_args)
    if pre_args.conf:
        for action in parser._actions:
            if action.dest in config_args:
                action.required = False

    return parser.parse_args()


def main():
    # parse the config file
    args = get_arguments()
    srr = args.srr
    global only_download
    only_download = args.only_download
    global skip_download
    skip_download = args.skip_download
    global remove_fastq
    remove_fastq = args.remove_fastq
    global remove_bam
    remove_bam = args.remove_bam
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
    if args.download:
        download_path = args.download
    else:
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
    # process srr
    process(srr)

if __name__ == "__main__":
    main()



