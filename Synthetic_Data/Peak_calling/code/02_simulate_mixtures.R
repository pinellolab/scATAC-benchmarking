library(data.table)
library(dplyr)

meta <- read.table("../data/simulation_1_design.tsv")
colnames(meta) <- c("simulation", "minor_population", "percent")

bcell <- fread("/data/aryee/caleb/hemeATAC/Bcell_reads.bed", sep = "\t", col.names = c("chr", "start", "end", "readid", "qual", "strand"), stringsAsFactors = TRUE)
cd4 <- fread("/data/aryee/caleb/hemeATAC/CD4_reads.bed",  sep = "\t", col.names = c("chr", "start", "end", "readid", "qual", "strand"), stringsAsFactors = TRUE)
mono <- fread("/data/aryee/caleb/hemeATAC/Mono_reads.bed",  sep = "\t", col.names = c("chr", "start", "end", "readid", "qual", "strand"), stringsAsFactors = TRUE)

n_reads <- 100000000 
n_bcell <- dim(bcell)[1]
n_cd4 <- dim(cd4)[1]
n_mono <- dim(mono)[1]

makeSyntheticBedsForMacs2 <- function(sim_id, minor_celltype, minor_percentage, seed = 10){
  
  set.seed(seed)
  n_reads_minor <- minor_percentage * 0.01 * n_reads
  n_reads_major <- (100-minor_percentage)/2 * 0.01 * n_reads
  
  n_reads_bcell <- n_reads_major
  n_reads_cd4 <- n_reads_major
  n_reads_mono <- n_reads_major
  
  # Establish minor population
  if(minor_celltype == "Bcell") n_reads_bcell <- n_reads_minor
  if(minor_celltype == "CD4") n_reads_cd4 <- n_reads_minor
  if(minor_celltype == "Mono") n_reads_mono <- n_reads_minor
  
  # Subsample
  df_bcell <- bcell[sample(1:n_bcell, n_reads_bcell),]
  df_cd4 <- cd4[sample(1:n_cd4, n_reads_cd4),]
  df_mono <- mono[sample(1:n_mono, n_reads_mono), ]
  
  # Assign minor population
  if(minor_celltype == "Bcell") minor_df <- df_bcell
  if(minor_celltype == "CD4") minor_df <- df_cd4
  if(minor_celltype == "Mono") minor_df <- df_mono
  
  write.table(rbind(df_bcell, df_cd4, df_mono), file = paste0(sim_id, "_",as.character(seed), "_mix.bed"),
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  
  write.table(minor_df, file = paste0(sim_id,"_",as.character(seed), "_minor.bed"),
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  sim_id
}

lapply(1:5, function(seed){
  lapply(1:15, function(i){
    makeSyntheticBedsForMacs2(as.character(meta[i,1]),
                              as.character(meta[i,2]),
                              as.numeric(meta[i,3]), seed = seed)
  })
})
