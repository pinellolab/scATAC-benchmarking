#!/bin/bash
OUTPUT='./sc-bams_barcodes'
mkdir -p $OUTPUT
FILEPATH='../../input/sc-bams_nodup/'
METAFILE='./SnapATAC_metadata.tsv'
{
	read
	while IFS=$'\t' read -r -a myArray
	do
		echo "${myArray[0]}"
		echo "${myArray[1]}"
		# echo $FILEPATH${myArray[1]}.st.bam
		# echo $OUTPUT/${myArray[1]}.with_barcode.bam
	./prependBarcode $FILEPATH${myArray[1]}.dedup.st.bam ${myArray[0]} $OUTPUT/${myArray[1]}.with_barcode.bam
	done
}< $METAFILE



