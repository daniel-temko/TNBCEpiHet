source("../utilities/utils.R")
source("../utilities/files.R")

annotate_heat_peaks <- function(prefix, 
                                n_groups, 
                                ann_df = prrx1_chip_peaks_ann, 
                                in_path = "chip_seq_prrx1_annotation/", 
                                out_path = "analysis/prrx1_chip_seq/summary/"){
  for(i in 1:n_groups){
    message(i)
    regions_f_name <- paste0(prefix, "_heat_peaks", i, ".bed")
    regions_ann_f_name <- paste0(prefix, "_heat_peaks", i, ".csv")
    regions_df <- read.table(file.path(in_path, regions_f_name))
    region_ids <- region_ids_from_zero_based_half_open(regions_df[,1], 
                                                       regions_df[,2], 
                                                       regions_df[,3])
    stopifnot(all(region_ids %in% rownames(ann_df)))
    regions_ann_df <- data.frame(region = pub_region_ids_from_region_ids(region_ids),
                                 gene = ann_df[region_ids, "gene.symbol"])
    write.table(regions_ann_df, file.path(out_path, regions_ann_f_name), 
                sep = ",", row.names = FALSE)
  }
}

setwd("../")

load("R_Data/PRRX1_ChIP_Data.RData")
load("R_Data/PRRX1allHeatmap_Data.RData")
all_n_groups <- n_groups
load("R_Data/PRRX1basHeatmap_Data.RData")
bas_n_groups <- n_groups
load("R_Data/PRRX1mesHeatmap_Data.RData")
mes_n_groups <- n_groups
load("R_Data/PRRX1neuroHeatmap_Data.RData")
neuro_n_groups <- n_groups
load("R_Data/PRRX1rhabdoidHeatmap_Data.RData")
rhabdoid_n_groups <- n_groups

annotate_heat_peaks("all", all_n_groups)
annotate_heat_peaks("bas", bas_n_groups)
annotate_heat_peaks("mes", mes_n_groups)
annotate_heat_peaks("neuro", neuro_n_groups)
annotate_heat_peaks("rhabdoid", rhabdoid_n_groups)