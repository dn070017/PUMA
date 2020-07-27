#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N read_and_print_networks
#$ -j y
#$ -o read_and_print_networks.qlog

echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "Task index number : $SGE_TASK_ID"
echo "=========================================================="

mkdir txt
matlab -nodisplay -nodesktop -nojvm -nosplash < read_and_print_networks.m > read_and_print_networks.qlog

echo "=========================================================="
echo "Finished on : $(date)"
echo "=========================================================="
