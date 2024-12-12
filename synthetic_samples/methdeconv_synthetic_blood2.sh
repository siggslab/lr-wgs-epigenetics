#!/bin/bash
#$ -S /bin/bash
#$ -o $HOME/jobs_logs
#$ -e $HOME/jobs_logs
#$ -cwd
#$ -q short.q
#$ -N syntehtic_methdeconv
#$ -l mem_requested=16G
#$ -pe smp 8

source ~/.bashrc && conda activate bamtools
cd /share/ScratchGeneral/yvefon/methylDeconv/ShreeData_synthetic_samples/configs

config=$(sed -n "${SGE_TASK_ID}p" ./configs_blood2.txt)

python ~/scripts/synthetic_sample_methatlas.py "${config}"
