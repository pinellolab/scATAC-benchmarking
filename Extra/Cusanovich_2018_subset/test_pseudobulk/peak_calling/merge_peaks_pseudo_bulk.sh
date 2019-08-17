############################################################
#Script to create pseudobulk peaks:                        #
# input: a file with column1 cellid and column2  clusterid #
#output: a peak file for every cluster                     #
#                                                          #
#Autor: Tommaso Andreani 08/01/2019                        #
############################################################


source activate biotools_huidong

#Format the file name
ls ../../input/sc-bams_nodup/ > files_name.txt
paste files_name.txt pseudobulk.tsv | cut -f 1,3 > pseudobulk_name.tsv 



echo 'start! assemble bam files for each cluster'
#Extract file name for each cluster
for i in $(seq 1 13);
do
	awk -F $'\t' -v i="$i" '$2 == i' pseudobulk.tsv | cut -f 1 > cluster.$i;
done

#Add path
for i in cluster* ;
do
  	perl -ne 'print "../../input/sc-bams_nodup/" . $_' $i > $i.input_path;
done

#Select only the file name
for i in *input_path;
do
	cut -f 1 $i > $i.file_name;
done


echo 'merge the bam files into one bam'
#Merge
for i in *.file_name; 
do 
	samtools merge $i.merged.bam -b $i ;
done


#Create index for the files
ls *bam | awk -F "." '{print $1"."$2".""input_path.file_name.merged"}' > bam_files.txt
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

