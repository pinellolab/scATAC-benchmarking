library(data.table)
library(BuenColors)
library(dplyr)
library(GenomicRanges)

"%ni%" <- Negate("%in%")

meta <- read.table("../data/simulation_1_design.tsv")

bcell_gr <- fread(paste0("zcat < ", "../peak_calls/Bcell_allData_peaks.narrowPeak.gz")) %>%
  data.frame() %>% makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3")
cd4_gr <- fread(paste0("zcat < ","../peak_calls/CD4_allData_peaks.narrowPeak.gz")) %>%
  data.frame() %>% makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3")
mono_gr <- fread(paste0("zcat < ","../peak_calls/Mono_allData_peaks.narrowPeak.gz")) %>%
  data.frame() %>% makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3")

ov_bc <- findOverlaps(bcell_gr, cd4_gr)
ov_bm <- findOverlaps(bcell_gr, mono_gr)
ov_cm <- findOverlaps(cd4_gr, mono_gr)

cts_b <- bcell_gr[1:length(bcell_gr) %ni% c(queryHits(ov_bc), queryHits(ov_bm))]
cts_c <- cd4_gr[1:length(cd4_gr) %ni% c(subjectHits(ov_bc), queryHits(ov_cm))]
cts_m <- mono_gr[1:length(mono_gr) %ni% c(subjectHits(ov_bm), subjectHits(ov_bm))]

eval_sim <- function(idx,sim){
  sim_id <- meta[idx,1] %>% as.character()
  minor_celltype <- meta[idx,2] %>% as.character()
  
  # Pull out the rate at which peaks are called from a clustered minor population
  minor_gr <- fread(paste0("zcat < ","../peak_calls/",sim_id, "_", as.character(sim),"_minor.bed_peaks.narrowPeak.gz")) %>%
    data.frame() %>% makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3")
  mix_gr <- fread(paste0("zcat < ","../peak_calls/",sim_id, "_", as.character(sim),"_mix.bed_peaks.narrowPeak.gz")) %>%
    data.frame() %>% makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3")
  ov <- findOverlaps(minor_gr, mix_gr)
  recovery_rate_all <- mean(1:length(minor_gr) %in% queryHits(ov))
  
  # Pull out cell-type-specific peaks
  if(minor_celltype == "Bcell") cts_gr <- cts_b
  if(minor_celltype == "CD4") cts_gr <- cts_c
  if(minor_celltype == "Mono") cts_gr <- cts_m
  
  ov_cts <- findOverlaps(cts_gr, minor_gr)
  cts_cc <- minor_gr[unique(subjectHits(ov_cts))]
  
  ov <- findOverlaps(cts_cc, mix_gr)
  recovery_rate_cts <- mean(1:length(cts_cc) %in% queryHits(ov))
  
  data.frame(meta[idx,], recovery_rate_all, recovery_rate_cts)
}


lapply(1:5, function(simn){
  lapply(1:15, eval_sim, sim = simn) %>% rbindlist() -> raw_df_onesim
  raw_df_onesim
}) %>% rbindlist() %>% data.frame() -> raw_df

raw_df %>%
  group_by(V2, V3) %>%
  summarize(mean = mean(recovery_rate_cts),
            sem = sd(recovery_rate_cts)) -> plot_df
#plot_df$sem  <- 0.1

threeb <- c("CD4" = "#0081C9", "Bcell" = "#BA7FD0", "Mono" = "#FF5A00")
pA  <- ggplot(plot_df, aes(x= V3, y =  mean*100, color = V2)) +
  geom_point()  + geom_line() +
  pretty_plot() + L_border() + 
  geom_errorbar(aes(ymin=mean*100-sem*100, ymax=mean*100+sem*100), width=.2)  +
  scale_color_manual(values  =  threeb) +
  labs(x = "Minor population prevalence (%)", y = "% cell type-specific peaks recovered",
       color = "Minor\npopulation")

cowplot::ggsave(pA, filename = "../output/final_plot_peaksim.pdf", width = 4, height = 2.7)
