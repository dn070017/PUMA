#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N combine_networks
#$ -j y
#$ -o combine_networks.qlog
echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "Task index number : $SGE_TASK_ID"
echo "=========================================================="

R --vanilla < combine_networks.R > combine_networks.Rout

echo "=========================================================="
echo "Finished on : $(date)"
echo "=========================================================="
