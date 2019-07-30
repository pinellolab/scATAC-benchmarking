#!/bin/bash

snaptools snap-pre  \
    --input-file=sciATAC_mouse.merged.sorted.bam  \
    --output-snap=cusanovich2018.snap  \
    --genome-name=mm9  \
    --genome-size=../../input/mm9/mm9.chrom.sizes  \
    --min-mapq=10  \
    --min-flen=0  \
    --max-flen=1000  \
    --keep-chrm=TRUE  \
    --keep-single=False  \
    --keep-secondary=False  \
    --overwrite=True  \
    --max-num=1000000  \
    --min-cov=0  \
    --verbose=True

snaptools snap-add-bmat    \
    --snap-file=cusanovich2018.snap \
    --bin-size-list 5000    \
    --verbose=True


