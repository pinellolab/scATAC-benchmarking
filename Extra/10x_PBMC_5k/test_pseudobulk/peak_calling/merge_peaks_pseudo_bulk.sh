############################################################
#Script to create pseudobulk peaks:                        #
# input: a file with column1 cellid and column2  clusterid #
#output: a peak file for every cluster                     #
#                                                          #
#Author: Tommaso Andreani 08/01/2019                        #
############################################################


source activate biotools_huidong

echo 'format the file'
cat pseudobulk.tsv | cut -f 1 | awk '$0=$0".dedup.st.bam"' > add_extension.txt;
cat add_extension.txt | perl -ne 'print "atac_v1_pbmc_5k_possorted_bam." . $_' $i  > add_extension_in_front.txt;
paste add_extension_in_front.txt pseudobulk.tsv | cut -f 1,3 > pseudobulk_name_extended.tsv;
rm add_extension*


echo 'start! assemble bam files for each cluster'
#Extract file name for each cluster
for i in $(seq 1 8);
do
	awk -F $'\t' -v i="$i" '$2 == i' pseudobulk_name_extended.tsv  > cluster.$i;
done

#Extract bam file name
for i in cluster*;
do
	cut -f 1 $i > $i.name;
done


#Add path extension
for i in *name ;
do
  	perl -ne 'print "../../input/sc-bams_nodup/" . $_' $i > $i.input_path;
done

echo 'merge the bam files into one bam'
#Merge
for i in *.input_path; 
do 
	samtools merge $i.merged.bam -b $i ;
done


#Create index for the files
ls *bam | awk -F "." '{print $1"."$2".""name.input_path.merged"}' > bam_files.txt
file_list=`cat bam_files.txt`

#Create folder where store the macs output and processes
mkdir macs;
cd macs/;


echo 'call peaks now...'
#Call peaks
for i in $file_list;
do
	macs2 callpeak -t ../$i.bam -f BAM -g 2.7e9 -n $i.merged --outdir $i --nomodel --shift -100 --extsize 200 2> $i.log;
done


#Clean
cd ../;
mkdir bam;
mv *bam bam/;

mkdir intermediate_files;
mv bam_files.txt samples* intermediate_files/;

echo "you are set the script is done ;) "

