""" 
Split a single bam file into single-cell bam files based on barcodes

@author Huidong Chen

"""

###  be sure to first sort bam file by cell barcode 'CB'  
### `samtools sort -t CB -@ 20 possorted_bam.dedup.bam > possorted_bam.dedup.sorted.bam`


import pandas as pd
import time
import pysam
import os

df_barcodes = pd.read_csv("./filtered_peak_bc_matrix/barcodes.tsv",header=None)

### Input varibles
# bam file to split
unsplit_file = "./possorted_bam.dedup.sorted.bam"
# directory of output files
out_dir = "./sc_bam_dedup/"


if(not os.path.exists(out_dir)):
        os.makedirs(out_dir)

# variable to hold barcode index
CB_hold = 'unset'
itr = 0
# read in upsplit file and process reads line by line
samfile = pysam.AlignmentFile(unsplit_file, "rb")

start_time = time.time()
for read in samfile.fetch(until_eof=True):
    if(len([ x for x in read.tags if x[0] == "CB"])>0):
        # barcode itr for current read
        CB_itr = read.get_tag( 'CB')
        if(CB_itr in df_barcodes[0].tolist()):
            # if change in barcode or first line; open new file  
            if(CB_itr!=CB_hold or itr==0):
                # close previous split file, only if not first read in file
                if(itr!=0):
                    split_file.close()
                CB_hold = CB_itr
                itr+=1
                split_file = pysam.AlignmentFile(out_dir + "{}.bam".format(CB_hold), "wb", template=samfile)
            # write reads with the same barcode to file
            split_file.write(read)    
split_file.close()
samfile.close()
end_time = time.time()


print(end_time - start_time)