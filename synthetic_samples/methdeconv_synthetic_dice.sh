#!/bin/bash
source ~/.bashrc && conda activate bamtools
# ^ conda env with bamtools

cd "<INSERT CONFIG FILE DIRECTORY>"

config=$(sed -n "${SGE_TASK_ID}p" ./configs.txt)
# ^ Text file, listing config files to process, with 1 file per line
# $SGE_TASK_ID is from the job submission system on the compute cluster used to run this
# it's basically an index, and can be replaced with an equivalent for loop over the lines in configs.txt

python synthetic_sample_methatlas.py "${config}"
