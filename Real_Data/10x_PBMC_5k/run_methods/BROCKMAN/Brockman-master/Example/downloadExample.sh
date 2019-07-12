#!/bin/bash

mkdir -p Downloaded_fastQs
cd Downloaded_fastQs
while read -r id
do
	fastq-dump --split-files --gzip  $id
done < ../SRA_entries.txt
