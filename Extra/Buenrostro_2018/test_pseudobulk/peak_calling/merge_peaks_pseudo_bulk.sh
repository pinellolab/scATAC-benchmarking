############################################################
#Script to create pseudobulk peaks:                        #
# input: a file with column1 cellid and column2  clusterid #
#output: a peak file for every cluster                     #
#                                                          #
#Autor: Tommaso Andreani 08/01/2019                        #
############################################################


source activate biotools_huidong

echo 'start! assemble bam files for each cluster'
#Extract file name for each cluster
for i in $(seq 1 10);
do
	awk -F $'\t' -v i="$i" '$2 == i' pseudobulk.tsv | cut -f 1 > cluster.$i;
done


echo 'add the right extension to the files'
#Add extension
for i in cluster.*;
do
	awk '$0=$0".dedup.st.bam"' $i > samples.$i;
done

rm cluster*


echo 'add the path to the bam'
#Add path
for i in samples.cluster.* ;
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
ls *bam | awk -F "." '{print $1"."$2"."$3}' > bam_files.txt
file_list=`cat bam_files.txt`

#Create folder where store the macs output and processes
mkdir macs;
cd macs/;


echo 'call peaks now...'
#Call peaks
for i in $file_list;
do
	macs2 callpeak -t ../$i.input_path.merged.bam -f BAM -g 2.7e9 -n $i.merged --outdir $i --nomodel --shift -100 --extsize 200 2> $i.log;
done


#Clean
cd ../;
mkdir bam;
mv *bam bam/;

mkdir intermediate_files;
mv bam_files.txt samples* intermediate_files/;

echo "you are set the script is done ;) "

