cd ~/github/MassiveQC/example
mkdir -p results/
mkdir -p results/pbs
conda activate MassiveQC

function ProcessOne {
  SingleQC -c example_conf.txt -s $1 --only_download
  # conda activate in pbs or slurm will be failed, see https://github.com/conda/conda/issues/5071
  echo "source activate MassiveQC; cd ~/github/MassiveQC/example; SingleQC -c example_conf.txt -s $1 --remove_fastq --remove_bam" > results/pbs/${1}.pbs
  qsub -q batch -V -l nodes=1:ppn=2 results/pbs/${1}.pbs -o results/pbs/ -e results/pbs/
}
export -f ProcessOne
sed "/^ *#/d" input.txt|sed "1d"|cut -f2 |parallel --ungroup --verbose -j 3 ProcessOne {}

# the feature data file in /home/mwshi/project/MassiveQC/Features
# the filter result file is ~/project/MassiveQC/result.csv
#
IsoDetect -i input.txt -o results/

source activate MassiveQC
cd ~/github/MassiveQC/example
MultiQC -a ~/miniconda3/etc/asperaweb_id_dsa.openssh \
 -i input.txt --skip_download \
 -f ~/project/FastQ_Screen_Genomes/fastq_screen.conf \
 -g dmel-all-r6.11.gtf \
 -x hisat-index/dmel.genome \
 -k splicesites.dmel.txt \
 -p picard.jar \
 -r dmel_r6-11.refflat \
 -w 4 -t 6 \
 -d results/download \
 -o result2




