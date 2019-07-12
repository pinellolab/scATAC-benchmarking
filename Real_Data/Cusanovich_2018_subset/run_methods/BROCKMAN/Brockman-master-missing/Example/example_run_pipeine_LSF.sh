#!/bin/bash

#BSUB -J BROCKMAN[1-366]
#BSUB -o run_lsf_pipe.log
#BSUB -e run_lsf_pipe.err
#BSUB -n 2
#BSUB -q short


echo Job $LSB_JOBID:$LSB_JOBINDEX started on $LSB_SUB_HOST: `date` #print to the output file the current job id and task, and host
source activate ATACseq_BROCKMAN3
source config.sh #load config file (just to get $logDir and $sampleFile)

# export id=`awk "NR==$LSB_JOBINDEX" $sampleFile` # input the $curIndexth line of this file into $id
# splitID=($id) #split id on whitespace into array that can be accessed like ${splitID[0]}
# export id=${splitID[0]} # First column of the sample file sheet; the sample ID (must be unique)
# export logPre=$logDir/$id
#redirect stdout and stderr to the log file
# exec 1>>$logPre.olog
# exec 2>&1
../brockman_pipeline config.sh $LSB_JOBINDEX

bjobs $LSB_JOBID | grep "^usage *$LSB_JOBINDEX:" #display job usage stats

echo running at $(date)
