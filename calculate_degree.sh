#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N calculate_degree
#$ -j y
#$ -o calculate_degree.qlog
echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "Task index number : $SGE_TASK_ID"
echo "=========================================================="

R --vanilla < calculate_degree.R > calculate_degree.Rout

echo "=========================================================="
echo "Finished on : $(date)"
echo "=========================================================="
