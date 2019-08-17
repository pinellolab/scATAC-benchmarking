#!/bin/bash

#BSUB -J count_reads[1-5335]
#BSUB -o count_reads_windows.out
#BSUB -e count_reads_windows.err
#BSUB -We 5
#BSUB -q vshort

source activate ATACseq_preprocess
bampath=../input/sc-bams_nodup/
dirlist=(`ls $bampath*.bam`)
# echo ${dirlist[$LSB_JOBINDEX-1]}
mkdir -p count_reads_windows_output
echo ./count_reads_windows_output/$(basename ${dirlist[$LSB_JOBINDEX-1]}).windows.txt
bedtools coverage -a hg19.windows.5kb.bed -b ${dirlist[$LSB_JOBINDEX-1]} > ./count_reads_windows_output/$(basename ${dirlist[$LSB_JOBINDEX-1]}).windows.txt