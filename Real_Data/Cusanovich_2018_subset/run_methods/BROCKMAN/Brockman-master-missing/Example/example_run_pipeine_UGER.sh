#!/bin/bash -l

echo Job $JOB_ID:$SGE_TASK_ID started on $HOST: `date` #print to the output file the current job id and task, and host
use Anaconda3
source activate BrockmanEnv
source config.sh #load config file (just to get $logDir and $sampleFile)

export id=`awk "NR==$SGE_TASK_ID" $sampleFile` # input the $curIndexth line of this file into $id
splitID=($id) #split id on whitespace into array that can be accessed like ${splitID[0]}
export id=${splitID[0]} # First column of the sample file sheet; the sample ID (must be unique)
export logPre=$logDir/$id
#redirect stdout and stderr to the log file
exec 1>>$logPre.olog
exec 2>&1
../brockman_pipeline config.sh $SGE_TASK_ID

qstat -j $JOB_ID | grep "^usage *$SGE_TASK_ID:" #display job usage stats
