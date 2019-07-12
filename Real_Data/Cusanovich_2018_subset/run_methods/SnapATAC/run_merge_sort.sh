#!/bin/bash

prefix='list_bamfiles_'
dirlist=(`ls ${prefix}*.txt`)

for FILE in ${dirlist[*]}
do
    echo $FILE
    ID=$(echo $FILE | awk -F "_" '{print $3}' | awk -F "." '{print $1}')
    samtools merge -@ 10 sciATAC_mouse.merged_${ID}.bam -b $FILE
    samtools sort -@ 10 -n sciATAC_mouse.merged_${ID}.bam -o sciATAC_mouse.merged.sorted_${ID}.bam
done

samtools merge -@ 10 sciATAC_mouse.merged.bam sciATAC_mouse.merged.sorted_*.bam
samtools sort -@ 10 -n sciATAC_mouse.merged.bam -o sciATAC_mouse.merged.sorted.bam


