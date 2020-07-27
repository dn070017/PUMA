#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N combine_and_degree
#$ -j y
#$ -o combine_and_degree.qlog
echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "Task index number : $SGE_TASK_ID"
echo "=========================================================="

# bash lioness_run.sh 1 -1
bash read_and_print_networks.sh
bash combine_networks.sh
bash calculate_degree.sh

echo "=========================================================="
echo "Finished on : $(date)"
echo "=========================================================="
