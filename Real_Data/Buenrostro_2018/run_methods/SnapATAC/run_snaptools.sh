#!/bin/bash

snaptools snap-pre  \
    --input-file=Buenrostro2018.merged.sorted.bam  \
    --output-snap=Buenrostro2018.snap  \
    --genome-name=hg19  \
    --genome-size=../../input/hg19/hg19.chrom.sizes  \
    --min-mapq=30  \
    --min-flen=0  \
    --max-flen=1000  \
    --keep-chrm=TRUE  \
    --keep-single=False  \
    --keep-secondary=False  \
    --overwrite=True  \
    --max-num=1000000  \
    --min-cov=100  \
    --verbose=True

snaptools snap-add-bmat    \
    --snap-file=Buenrostro2018.snap \
    --bin-size-list 5000    \
    --verbose=True


