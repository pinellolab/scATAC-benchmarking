# scATAC-benchmarking

Recent innovations in single-cell Assay for Transposase Accessible Chromatin using sequencing (scATAC-seq) enable profiling of the epigenetic landscape of thousands of individual cells. scATAC-seq data analysis presents unique methodological challenges. scATAC-seq experiments sample DNA, which, due to low copy numbers (diploid in humans) lead to inherent data sparsity (1-10% of peaks detected per cell) compared to transcriptomic (scRNA-seq) data (20-50% of expressed genes detected per cell). Such challenges in data generation emphasize the need for informative features to assess cell heterogeneity at the chromatin level.  

<img src="images/Figure1.png">

We present a benchmarking framework that was applied to 10 computational methods for scATAC-seq on 13 synthetic and real datasets from different assays, profiling cell types from diverse tissues and organisms. Methods for processing and featurizing scATAC-seq data were evaluated by their ability to discriminate cell types when combined with common unsupervised clustering approaches. We rank evaluated methods and discuss computational challenges associated with scATAC-seq analysis including inherently sparse data, determination of features, peak calling, the effects of sequencing coverage and noise, and clustering performance. Running times and memory requirements are also discussed. 

<img src="images/Figure2.png">


Please check out our preprint in bioRxiv: [Chen, H. et al. Assessment of computational methods for the analysis of single-cell ATAC-seq data]()

Single Cell ATAC-seq Benchmarking Framework
-------------------------------------------

* Real Data

- The _**Buenrostro2018**_ dataset (folder name **'Buenrostro_2018'**)

  - input
  - run_methods
    - [BROCKMAN](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/run_methods/BROCKMAN/BROCKMAN_buenrostro2018.ipynb?flush_cache=true)
    - [chromVAR-kmers](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/run_methods/chromVAR/chromVAR_buenrostro2018_kmers.ipynb?flush_cache=true)
    - [chromVAR-motifs](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/run_methods/chromVAR/chromVAR_buenrostro2018_motifs.ipynb?flush_cache=true)
    - [Cicero](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/run_methods/Cicero/Cicero_buenrostro2018.ipynb?flush_cache=true)
    - [cisTopic](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/run_methods/Cicero/Cicero_buenrostro2018.ipynb?flush_cache=true)
    - [Control](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/run_methods/Control/Control_buenrostro2018.ipynb?flush_cache=true)
    - [Cusanovich2018](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/run_methods/Cusanovich2018/Cusanovich2018_buenrostro2018.ipynb?flush_cache=true)
    - [GeneScoring](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/run_methods/GeneScoring/GeneScoring_buenrostro2018.ipynb?flush_cache=true)
    - [scABC](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/run_methods/scABC/scABC_buenrostro2018.ipynb?flush_cache=true)
    - [Scasat](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/run_methods/Scasat/Scasat_buenrostro2018.ipynb?flush_cache=true)
    - [SCRAT](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/run_methods/SCRAT/SCRAT_buenrostro2018.ipynb?flush_cache=true)
    - [SnapATAC](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/run_methods/SnapATAC/SnapATAC_buenrostro2018.ipynb?flush_cache=true)
  - [clustering](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/run_clustering_buenrostro2018.ipynb?flush_cache=true)
  - [UMAP visualization](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/run_umap_buenrostro2018.ipynb?flush_cache=true)
  - end-to-end clustering (folder name 'extra_clustering')
    - [Cicero](hhttps://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/extra_clustering/Cicero/Cicero_buenrostro2018.ipynb?flush_cache=true)
    - [cisTopic](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/extra_clustering/cisTopic/cisTopic_buenrostro2018.ipynb?flush_cache=true)
    - [Cusanovich2018](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/extra_clustering/Cusanovich2018/Cusanovich2018_buenrostro2018.ipynb?flush_cache=true)
    - [scABC](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/extra_clustering/scABC/scABC_buenrostro2018.ipynb?flush_cache=true)
    - [Scasat](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/extra_clustering/Scasat/Scasat_buenrostro2018.ipynb?flush_cache=true)
    - [SnapATAC](https://nbviewer.jupyter.org/github/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/extra_clustering/SnapATAC/SnapATAC_buenrostro2018.ipynb?flush_cache=true)

- The _**Buenrostro2018 using bulk peaks**_ dataset (folder name **'Buenrostro_2018_bulkpeaks'**)


- The _**10X PBMCs**_ dataset (folder name **'10x_PBMC_5k'**)


- The _**downsampled sci-ATAC-seq-mouse**_ dataset (folder name **'Cusanovich_2018_subset'**)


- The _**full sci-ATAC-seq-mouse**_ dataset (folder name **'Cusanovich_2018'**)



* Synthetic Data



* Extra


