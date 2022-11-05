# MassiveQC

## Overview

The MassiveQC is a python-based tool for quality control of massive RNA-seq data. It mainly consists of the following parts
* Interpret input files and download sequencing files
* Check fastq files, filter unqualified read and fastq files.
* Align all sample batches and use multiple tools to analyze sample features.
* Identify outliers using isolation forests on sample features.

## Dependencies
MassiveQC is implemented in [Python 3.10](https://www.python.org/downloads/).


### System requirements
* Linux or Mac  
* Python 3.7 or later

### Prerequisites
MassiveQC has the following dependencies:

Python package:
* [scikit-learn](https://scikit-learn.org/), a Python module for machine learning.
* [shap](https://github.com/slundberg/shap), a game theoretic approach to explain the output of any machine learning model.
* [pysradb](https://github.com/saketkc/pysradb), a Python package for retrieving metadata and downloading datasets from SRA/ENA/GEO.
* [xopen](https://__github__.com/pycompression/xopen), a Python-based package for fast processing of compressed files..  
* [NumPy](https://numpy.org/), a fundamental package for scientific computing with Python.    
* [Pandas](https://pandas.pydata.org/) >= 1.3.5, a fast, powerful and flexible open source data analysis and manipulation tool.
* [fastparquet](https://github.com/dask/fastparquet/), a python implementation of the parquet format, aiming integrate into python-based big data work-flows.
* [more-itertools](https://pypi.org/project/more-itertools/), provides additional building blocks, recipes, and routines for working with Python iterables.
* [tqdm](https://github.com/tqdm/tqdm), a Fast, Extensible Progress Bar for Python and CLI.

Software:
* [Atropos](https://github.com/jdidion/atropos), atropos is tool for specific, sensitive, and speedy trimming of NGS reads.
* [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/), A tool for multi-genome mapping and quality control.
* [HISAT2](http://daehwankimlab.github.io/hisat2/), a fast and sensitive alignment program for mapping next-generation sequencing reads.
* [samtools](https://subread.sourceforge.net/featureCounts.html), a highly efficient general-purpose read summarization program.
* [bamtools](https://subread.sourceforge.net/featureCounts.html), a highly efficient general-purpose read summarization program.
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
To create a new environment based on the environment.yaml:
```
conda env create -f environment.yml
```
Or for a basic environment and downloading optional dependencies as needed:
```
conda create -n MassiveQC -c bioconda python=3 MassiveQC
```

## Usage
For users running locally, we provide the `MultiQC` to automate parallel processing.
```
usage: MultiQC [-h] [-c CONF] -i INPUT [-a ASCP_KEY] -f FASTQ_SCREEN_CONFIG -g GTF -x HT2_IDX [-k KNOWN_SPLICESITE_INFILE]
               -p PICARD -r REF_FLAT -o OUTDIR [-w WORKERS] [-t THREADS] [-d DOWNLOAD] [--only_download] [--skip_download]
               [--remove_fastq] [--remove_bam]

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
  --skip_download       Skip the download step
  --remove_fastq        Don't remain the fastq after running hisat2
  --remove_bam          Don't remain the bam after running FeatureCounts
```

In the example, Users need to provide multiple files:
* `asperaweb_id_dsa.openssh` is the aspera key in [IBM aspera](https://www.ibm.com/products/aspera).
* `fastq_screen.conf` is the reference for [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/). It can be downloaded with `fastq_screen --get_genomes`.
* `dmel.genome` is Hisat2 index. Need to be created by the users.
* `dmel_r6-11.refflat` is created from gtf, Generated using the code below
* `picard.jar` can be download [here](https://github.com/broadinstitute/picard/releases)
```
# Generate refflat file
gtfToGenePred -genePredExt dmel-all-r6.11.gtf \
-ignoreGroupsWithoutExons /dev/stdout| awk 'BEGIN { OFS="\t"} \
{print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' >  dmel_r6-11.refflat
```
Users can optionally specify parameters using a conf file like `example/example_conf.txt`.

```
# Use conf file
MultiQC -c example/example_conf.txt -i example/input.txt

# Without conf file
MultiQC -a ~/miniconda3/etc/asperaweb_id_dsa.openssh \
-f ~/FastQ_Screen_Genomes/fastq_screen.conf \
-g example/dmel-all-r6.11.gtf \
-x example/dmel.genome \
-k example/splicesites.dmel.txt \
-p example/picard.jar \
-r example/dmel_r6-11.refflat \
-o ~/project/MassiveQC \
-i example/input.txt
```

For users running on cluster (PBS or Slurm), we provide `SingleQC` and `IsoDetect` for single sample. Users can batch process samples according to their platform.

```
usage: SingleQC [-h] [-c CONF] -s SRR [-a ASCP_KEY] -f FASTQ_SCREEN_CONFIG -g GTF -x HT2_IDX [-k KNOWN_SPLICESITE_INFILE] -p
                PICARD -r REF_FLAT -o OUTDIR [-t THREADS] [-d DOWNLOAD] [--only_download] [--skip_download] [--remove_fastq]
                [--remove_bam]

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
                        Path to SRA fastq files. The default is $OUTDIR/download
  --only_download       Only run the download step
  --skip_download       Skip the download step
  --remove_fastq        Don't remain the fastq after running hisat2
  --remove_bam          Don't remain the bam after running FeatureCounts
```

Here we provide an example on a PBS.
```
# create the outdir and pbs directory
cd ~/github/MassiveQC/example # my work directory
mkdir -p results/
mkdir -p results/pbs

# Generate pbs script and submit using for loop 
function ProcessOne {
  SingleQC -c example_conf.txt -s $1 --only_download
  # conda activate in pbs or slurm will be failed, see https://github.com/conda/conda/issues/5071
  echo "cd ~/github/MassiveQC/example; SingleQC -c example_conf.txt -s $1 --remove_fastq --remove_bam" > results/pbs/${1}.pbs
  qsub -q batch -V -l nodes=1:ppn=2 results/pbs/${1}.pbs -o results/pbs/ -e results/pbs/
}
export -f ProcessOne
sed "/^ *#/d" input.txt|sed "1d"|cut -f2 |parallel --ungroup --verbose -j 3 ProcessOne {}



# After all scripts have been run, start building the feature matrix
# Get started with outlier filtering
# the feature data file in results/Features
IsoDetect -i input.txt -o results/
# the filter result file is results/result.csv
```
We provide the test file in example directory
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
MassiveQC will process and feature extraction for each srr, then merge the features through srx, and finally identify qualified srx samples.

Users can also provide only srr, so that massiveQC will filter only on srr
```
srr
SRR6826933
SRR21851295
SRR21851296
SRR039433
```

### Output
MassiveQC will generate multiple result files under the OUTDIR folder, the path tree is as follows.
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

