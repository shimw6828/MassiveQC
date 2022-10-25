mkdir -p ~/project/MassiveQC
mkdir -p ~/project/MassiveQC/pbs
cd ~/project/MassiveQC/pbs
for srr in `sed "1d" ~/project/input.txt|cut -f2 |xargs`
do
  #conda activate in pbs or slurm will be failed, see https://github.com/conda/conda/issues/5071
  echo "source activate MassiveQC; SingleQC -c ~/project/config.txt -s $srr --skip_download" > ~/project/MassiveQC/pbs/${srr}.pbs
  qsub -q batch -V -l nodes=1:ppn=2 ~/project/MassiveQC/pbs/${srr}.pbs
done

IsoDetect -i ~/project/input.txt -o ~/project/MassiveQC/
#the feature data file in /home/mwshi/project/MassiveQC/Features
#the filter result file is ~/project/MassiveQC/result.csv



