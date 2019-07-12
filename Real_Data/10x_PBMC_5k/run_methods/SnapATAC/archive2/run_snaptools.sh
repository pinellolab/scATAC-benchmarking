#!/bin/bash

snaptools snap-pre  \
    --input-file=10xpbmc5k.snap.nsrt.bam  \
    --output-snap=10xpbmc5k.snap  \
    --genome-name=hg19  \
    --genome-size=../../input/hg19/hg19.chrom.sizes  \
    --min-mapq=10  \
    --min-flen=0  \
    --max-flen=1000  \
    --keep-chrm=TRUE  \
    --keep-single=False  \
    --keep-secondary=False  \
    --overwrite=True  \
    --max-num=20000  \
    --min-cov=0  \
    --verbose=True

snaptools snap-add-bmat    \
    --snap-file=10xpbmc5k.snap \
    --bin-size-list 5000    \
    --verbose=True