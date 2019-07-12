#!/bin/bash

#BSUB -J count_reads[1-2034]
#BSUB -o count_reads_peaks.out
#BSUB -e count_reads_peaks.err
#BSUB -We 5
#BSUB -q vshort

source activate ATACseq_preprocess
bampath=../../input/sc-bams_nodup/
dirlist=(`ls $bampath*.bam`)
# echo ${dirlist[$LSB_JOBINDEX-1]}
mkdir -p count_reads_peaks_output
echo ./count_reads_peaks_output/$(basename ${dirlist[$LSB_JOBINDEX-1]}).peaks.txt
bedtools coverage -a ../../input/GSE96769_PeakFile_20160207.bed -b ${dirlist[$LSB_JOBINDEX-1]} | cut -f 1,2,3,15- > ./count_reads_peaks_output/$(basename ${dirlist[$LSB_JOBINDEX-1]}).peaks.txt