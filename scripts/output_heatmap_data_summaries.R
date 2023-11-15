source("../utilities/preprocessing.R")
source("../utilities/utils.R")

setwd("../")

load(file = "R_Data/ChIP_Data.RData")
load(file = "R_Data/RNA_Data.RData")
load(file = "R_Data/Methyl_Data.RData")
load(file = "R_Data/Metab_Data.RData")
load(file = "R_Data/HMS_Data.RData")
load(file = "R_Data/BH3_Data.RData")
load(file = "R_Data/DS_Data.RData")

perc_features <- 20

message("Filtering")

chip_se_data <- get_heatmap_data(chip_se_all_rpkm_norm, perc_features = 20)
chip_se_data$gene <- chip_se_ann[rownames(chip_se_data), "gene.symbol"]
chip_se_data$superenhancer <- pub_region_ids_from_region_ids(rownames(chip_se_data))
chip_se_data <- bring_multiple_to_front(chip_se_data, c("superenhancer", "gene"))

rna_data <- get_heatmap_data(rna_rpkm_norm, perc_features = 20)
rna_data$gene <- rownames(rna_data)
rna_data <- bring_to_front(rna_data, "gene")

methyl_se_data <- get_heatmap_data(methyl_se_agg_m, perc_features = 20)
methyl_se_data$gene <- methyl_se_agg_ann[rownames(methyl_se_data), "gene"]
methyl_se_data$superenhancer <- pub_region_ids_from_region_ids(rownames(methyl_se_data))
methyl_se_data <- bring_multiple_to_front(methyl_se_data, c("superenhancer", "gene"))

hms_data <- get_heatmap_data(hms_norm, perc_features = 100)
hms_data$feature <- rownames(hms_data)
hms_data <- bring_to_front(hms_data, "feature")

metab_data <- get_heatmap_data(metab_norm, perc_features = 100)
metab_data$metabolite <- rownames(metab_data)
metab_data <- bring_to_front(metab_data, "metabolite")

ds_data <- get_heatmap_data(ds_auc, perc_features = 100)
ds_data$agent <- rownames(ds_data)
ds_data <- bring_to_front(ds_data, "agent")

bh3_data <- get_heatmap_data(bh3_conc, perc_features = 100)
bh3_data$feature <- rownames(bh3_data)
bh3_data <- bring_to_front(bh3_data, "feature")

message("Writing to file")

write.table(chip_se_data, 
            "analysis/chip_seq/summary/chip_seq_se_heatmap_data.csv", 
            sep = ",", 
            row.names = FALSE)
write.table(rna_data, 
            "analysis/rna_seq/summary/rna_seq_heatmap_data.csv", 
            sep = ",", 
            row.names = FALSE)
write.table(methyl_se_data, 
            "analysis/dna_methylation/summary/dna_methylation_se_heatmap_data.csv", 
            sep = ",", 
            row.names = FALSE)
write.table(hms_data, 
            "analysis/hms/summary/hms_heatmap_data.csv", 
            sep = ",", 
            row.names = FALSE)
write.table(metab_data, 
            "analysis/metabolomics/summary/metabolomics_heatmap_data.csv", 
            sep = ",", 
            row.names = FALSE)
write.table(ds_data,
            "analysis/drug_screen/summary/drug_screen_heatmap_data.csv", 
            sep = ",", 
            row.names = FALSE)
write.table(bh3_data, 
            "analysis/bh3_profiling/summary/bh3_profiling_heatmap_data.csv",
            sep = ",", 
            row.names = FALSE)
