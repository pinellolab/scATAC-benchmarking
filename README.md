# scATAC-benchmarking

Recent innovations in single-cell Assay for Transposase Accessible Chromatin using sequencing (scATAC-seq) enable profiling of the epigenetic landscape of thousands of individual cells. scATAC-seq data analysis presents unique methodological challenges. scATAC-seq experiments sample DNA, which, due to low copy numbers (diploid in humans) lead to inherent data sparsity (1-10% of peaks detected per cell) compared to transcriptomic (scRNA-seq) data (20-50% of expressed genes detected per cell). Such challenges in data generation emphasize the need for informative features to assess cell heterogeneity at the chromatin level.  

<img src="https://github.com/pinellolab/scATAC-benchmarking/images/Figure1.png">

We present a benchmarking framework that was applied to 10 computational methods for scATAC-seq on 13 synthetic and real datasets from different assays, profiling cell types from diverse tissues and organisms. Methods for processing and featurizing scATAC-seq data were evaluated by their ability to discriminate cell types when combined with common unsupervised clustering approaches. We rank evaluated methods and discuss computational challenges associated with scATAC-seq analysis including inherently sparse data, determination of features, peak calling, the effects of sequencing coverage and noise, and clustering performance. Running times and memory requirements are also discussed. 

<img src="https://github.com/pinellolab/scATAC-benchmarking/images/Figure2.png">


Please check out our preprint in BioRxiv: 

Single Cell ATAC-seq Benchmarking Framework
-------------------------------------------


