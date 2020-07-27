#!/bin/bash
# Working script for running PUMA in different scenarios.
# Edit 'puma_config.m' to set program parameters.

# Console running
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('puma_config.m'); run('puma_run.m');"

# Foreground running
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('puma_config.m'); run('puma_run.m'); quit;"

# Background running (./puma_run.sh &)
nohup matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('puma_config.m'); run('puma_run.m'); quit;" >& puma.`hostname`.log

# Email notification when done
echo "PUMA run on `hostname` has just finished: `date`." | mail -s "Task finished on `hostname`" -a puma.`hostname`.log `whoami`

# Testing
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('puma_config_test.m'); run('puma_run.m'); quit;"
#diff /tmp/puma.test.txt test_data/puma.txt

# Benchmark (5 repeats)
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('puma_config_test.m'); run('puma_run.m'); run('puma_run.m'); run('puma_run.m'); run('puma_run.m'); run('puma_run.m'); quit;" >& puma.test.`hostname`.log
