# MassiveQC

## Overview

The MassiveQC is a python-based tool for quality control of massive RNA-seq data. It mainly consists of the following parts
* Interpret input files and download sequencing files
* Check fastq files, filter unqualified read and fastq files.
* Align all sample batches and use multiple tools to analyze sample features.
* Identify Outliers Using Isolation Forests on Sample Features.

## Dependencies
MassiveQC is implemented in [Python 3.10](https://www.python.org/downloads/).


### System requirements
* Linux or Mac  
* Python 3.7 or later

### Prerequisites
MntJULiP has the following dependencies:
* [sklearn](https://scikit-learn.org/), a Python module for machine learning.
* [shap](https://github.com/slundberg/shap), a game theoretic approach to explain the output of any machine learning model.
* [pysradb](https://github.com/saketkc/pysradb), a Python package for retrieving metadata and downloading datasets from SRA/ENA/GEO.
* [xopen](https://github.com/pycompression/xopen), a Python-based package for fast processing of compressed files..  
* [NumPy](https://numpy.org/), a fundamental package for scientific computing with Python.    
* [Pandas](https://pandas.pydata.org/), a fast, powerful and flexible open source data analysis and manipulation tool.  
* [matplotlib](https://matplotlib.org/), a comprehensive library for creating static, animated, and interactive visualizations in Python.
* [seaborn](https://seaborn.pydata.org/), a Python data visualization library based on matplotlib.
* [Atropos](https://github.com/jdidion/atropos), atropos is tool for specific, sensitive, and speedy trimming of NGS reads.
* [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/), A tool for multi-genome mapping and quality control.
* [HISAT2](http://daehwankimlab.github.io/hisat2/), a fast and sensitive alignment program for mapping next-generation sequencing reads.
* [featureCounts](https://subread.sourceforge.net/featureCounts.html), a highly efficient general-purpose read summarization program.

To install using pip
```
pip3 install MassiveQC
```
Alternatively, if you use conda:
```
conda install -c bioconda MassiveQC
```
This step will install all the dependencies. If you have an existing environment with a lot of pre-installed packages, conda might be slow. Please consider creating a new enviroment for MassiveQC:
```
conda create -c bioconda -n MassiveQC MassiveQC
```

## Usage
For users running locally, we provide the `MultiQC` to automate parallel processing.
```
usage: MultiQC [-h] [-c CONF] -i INPUT [-a ASCP_KEY] -f FASTQ_SCREEN_CONFIG -g GTF -x HT2_IDX [-k KNOWN_SPLICESITE_INFILE] -p
                       PICARD -r REF_FLAT -o OUTDIR [-w WORKERS] [-t THREADS] [-d DOWNLOAD] [--only_download]

...

options:
  -h, --help            show this help message and exit
  -c CONF, --conf CONF
  -i INPUT, --input INPUT
                        Input file, containing two columns srx and srr
  -a ASCP_KEY, --ascp_key ASCP_KEY
                        Locate aspera key. Default $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh
  -f FASTQ_SCREEN_CONFIG, --fastq_screen_config FASTQ_SCREEN_CONFIG
                        Path to the fastq_screen conf file, can be download from fastq_screen website
  -g GTF, --gtf GTF     Path to the GTF file with annotations
  -x HT2_IDX, --ht2-idx HT2_IDX
                        Hisat2 index filename prefix
  -k KNOWN_SPLICESITE_INFILE, --known-splicesite-infile KNOWN_SPLICESITE_INFILE
                        Hisat2 splicesite file, provide a list of known splice sites
  -p PICARD, --picard PICARD
                        Path to picard.jar
  -r REF_FLAT, --ref_flat REF_FLAT
                        Path to refflat file
  -o OUTDIR, --outdir OUTDIR
                        Path to result output directory. If it doesn't exist, it will be created automatically
  -w WORKERS, --workers WORKERS
                        The number of simultaneous tasks
  -t THREADS, --THREADS THREADS
                        The number of threads for tools like Hisat2 in one task
  -d DOWNLOAD, --download DOWNLOAD
                        Path to SRA fastq files. The default is $OUTDIR/download
  --only_download       Only run the download step
```
Users can optionally specify parameters using a conf file like `example/example_conf.txt`.
`picard.jar` can be download [here](https://github.com/broadinstitute/picard/releases)
```
# Use conf file
MultiQC -c example/example_conf.txt -i example/srx2srr.txt

# Without conf file
MultiQC -a ~/miniconda3/etc/asperaweb_id_dsa.openssh \
-f ~/FastQ_Screen_Genomes/fastq_screen.conf \
-g example/dmel-all-r6.11.gtf \
-x example/dmel.genome \
-k example/splicesites.dmel.txt \
-p example/picard.jar \
-r example/dmel_r6-11.refflat \
-o ~/project/MassiveQC \
-i example/srx2srr.txt
```
In the example, Users need to provide multiple files:
* `asperaweb_id_dsa.openssh` is the aspera key in [IBM aspera](https://www.ibm.com/products/aspera).
* `fastq_screen.conf` is the reference for [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/). It can be downloaded with `fastq_screen --get_genomes`.
* `dmel.genome` is Hisat2 index. Need to be created by the users.
* `dmel_r6-11.refflat` is created from gtf. Generated using the code below
```
# Generate refflat file
gtfToGenePred -genePredExt dmel-all-r6.11.gtf \
-ignoreGroupsWithoutExons /dev/stdout| awk 'BEGIN { OFS="\t"} \
{print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' >  dmel_r6-11.refflat
```

For users running on cluster (PBS or Slurm), we provide `SingleQC` and `IsoDetect` for single sample. Users can batch process samples according to their platform.

```
usage: SingleQC [-h] [-c CONF] -s SRR [-a ASCP_KEY] -f FASTQ_SCREEN_CONFIG -g GTF -x
                        HT2_IDX [-k KNOWN_SPLICESITE_INFILE] -p PICARD -r REF_FLAT -o OUTDIR
                        [-t THREADS] [-d DOWNLOAD] [--only_download]

...

options:
  -h, --help            show this help message and exit
  -c CONF, --conf CONF
  -s SRR, --srr SRR     SRR id
  -a ASCP_KEY, --ascp_key ASCP_KEY
                        Locate aspera key. Default $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh
  -f FASTQ_SCREEN_CONFIG, --fastq_screen_config FASTQ_SCREEN_CONFIG
                        Path to the fastq_screen conf file, can be download from fastq_screen website
  -g GTF, --gtf GTF     Path to the GTF file with annotations
  -x HT2_IDX, --ht2-idx HT2_IDX
                        Hisat2 index filename prefix
  -k KNOWN_SPLICESITE_INFILE, --known-splicesite-infile KNOWN_SPLICESITE_INFILE
                        Hisat2 splicesite file, provide a list of known splice sites
  -p PICARD, --picard PICARD
                        Path to picard.jar
  -r REF_FLAT, --ref_flat REF_FLAT
                        Path to refflat file
  -o OUTDIR, --outdir OUTDIR
                        Path to result output directory. If it doesn't exist, it will be created automatically
  -t THREADS, --THREADS THREADS
                        The number of threads for tools like Hisat2 in one task
  -d DOWNLOAD, --download DOWNLOAD
                        Path to SRA fastq files. The default is $outdir/download
  --only_download       Only run the download step

```

Here we provide an example on a PBS.
```
# create the outdir and pbs dir
mkdir -p ~/project/MassiveQC
mkdir -p ~/project/MassiveQC/pbs
cd ~/project/MassiveQC/pbs
# Generate pbs script and submit using for loop 
for srr in `sed "1d" ~/project/input.txt|cut -f2 |xargs`
do
  # We have no network in the node, so we need to download the fastq file in advance
  SingleQC -c ~/project/config.txt -s $srr --only_download
  #conda activate in pbs or slurm will be failed, see https://github.com/conda/conda/issues/5071
  # If the download is successful, the download step will be skipped.
  echo "source activate MassiveQC; SingleQC -c ~/project/config.txt -s $srr --skip_download" > ~/project/MassiveQC/pbs/${srr}.pbs
  qsub -q batch -V -l nodes=1:ppn=2 ~/project/MassiveQC/pbs/${srr}.pbs
done

# After all scripts have been run, start building the feature matrix
# the result feature data file is in /home/mwshi/project/MassiveQC/Features
# Get started with outlier filtering
IsoDetect -f ~/project/MassiveQC/Features/features.parquet
#the filter result file is ~/project/MassiveQC/result.csv
```

## Input/Output
### Input
The main input of MassiveQC is a file contains a list of SRX and SRR id. For example:
```
srx srr
SRX3783112  SRR6826933
SRX17839866 SRR21851295
SRX17839866 SRR21851296
SRX018864   SRR039433
```
Massiveqc will process and feature extraction for each srr, then merge the features through srx, and finally identify qualified srx samples.

Users can also provide only srr, so that massiveQC will filter only on srr
```
srr
SRR6826933
SRR21851295
SRR21851296
SRR039433
```

### Output
MassiveQC will generate multiple result files under the outdir folder, the path tree is as follows.
```
OUTDIR
-Bam
-Count
-download # the downloaded fastq file
-QC_dir # the fastq file after quality control
-Feature # Sample features
 -aln_stats
 -atropos
 -count_summary
 -DoneSample
 -fastq_screen
 -genebody_coverage
 -hisat2
 -layout
 -markduplicates
 -rnaseqmetrics
 -strand
-result.csv # The result file, containing inlier and outlier samples.
```

