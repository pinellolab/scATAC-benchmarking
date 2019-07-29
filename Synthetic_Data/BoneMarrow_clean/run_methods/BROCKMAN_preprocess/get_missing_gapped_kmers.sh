#!/bin/bash

#BSUB -J brockman[1-1]
#BSUB -o brockman_missing.out
#BSUB -e brockman_missing.err
#BSUB -We 5
#BSUB -q vshort

source activate ATACseq_BROCKMAN3

dirlist=(`ls -d -1 ${PWD}/bed_missing/*`)
FILE=${dirlist[$LSB_JOBINDEX-1]}
id=$(basename "$FILE" | cut -d. -f1)
echo $id
mkdir -p bed_valid
bedtools intersect -a $FILE -b chroms.bed > ./bed_valid/$id.valid.bed
mkdir -p seqs
bedtools merge -i ./bed_valid/$id.valid.bed | awk 'BEGIN{FS="\t";OFS="\t"};{print $1":"$2"-"$3}' > ./seqs/$id.42bit
twoBitToFa ../../input/hg19/hg19.2bit  ./seqs/$id.fa -seqList=./seqs/$id.42bit
gzip ./seqs/$id.fa
mkdir -p kmers
AMUSED -bc -nsz -ds -ns -s 8 -q ./seqs/$id.fa.gz -o ./kmers/$id.scan
cat ./kmers/$id.scan | awk 'BEGIN{OFS="\t"; OFS="\t"};{print $1,$3}' | gzip -c > ./kmers/$id.freq.gz
