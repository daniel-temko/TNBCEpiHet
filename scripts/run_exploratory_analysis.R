library(ggplot2)
library(ggrepel)
library(umap)
library(pheatmap)

setwd("../")

source("scripts/utilities/exploratory_plots.R")
source("scripts/utilities/preprocessing.R")

# Load  data ----------------------------------------------------------------------------

load(file = "R_Data/Heatmap_Metadata.RData")

load(file = "R_Data/ChIP_Data.RData")
load(file = "R_Data/RNA_Data.RData")
load(file = "R_Data/Methyl_Data.RData")
load(file = "R_Data/Metab_Data.RData")
load(file = "R_Data/HMS_Data.RData")
load(file = "R_Data/BH3_Data.RData")
load(file = "R_Data/DS_Data.RData")

#exclusion_list <- list(hms = "HDQP1")

# Analyse RNA-seq data ----------------------------------------------------------------------------

rna_ranks <- variable_features(rna_rpkm_norm)

GenerateExploratoryPlots(data = rna_rpkm_norm,
                         top_features = rna_ranks,
                         n_pcs = rep(10, 6),
                         base_col = "#FF00FF",
                         base_dir = "analysis/rna_seq/exploratory/")

# Analyse ChIP-seq data ----------------------------------------------------------------------------

# Superenhancers
chip_se_ranks <- variable_features(chip_se_all_rpkm_norm)

GenerateExploratoryPlots(data = chip_se_all_rpkm_norm,
                         top_features = chip_se_ranks,
                         n_pcs = rep(10, 6),
                         base_col = "#FF8000",
                         base_dir = "analysis/chip_seq/exploratory/superenhancers/")

# Enhancers
chip_enh_all_ranks <- variable_features(chip_enh_all_rpkm_norm)

GenerateExploratoryPlots(data = chip_enh_all_rpkm_norm,
                         top_features = chip_enh_all_ranks,
                         n_pcs = rep(10, 6),
                         base_col = "#FF8000",
                         base_dir = "analysis/chip_seq/exploratory/enhancers/")

# Peaks
chip_peak_all_ranks <- variable_features(chip_peak_all_rpkm_norm)

GenerateExploratoryPlots(data = chip_peak_all_rpkm_norm,
                         top_features = chip_peak_all_ranks,
                         n_pcs = rep(10, 6),
                         base_col = "#FF8000",
                         base_dir = "analysis/chip_seq/exploratory/peaks/")


# Analyse DNA methylation data ----------------------------------------------------------------------------

# Intergenic super-enhancers
methyl_se_int_variances <- variable_features(methyl_se_int_agg_m)

GenerateExploratoryPlots(data = methyl_se_int_agg_m,
                         top_features = methyl_se_int_variances,
                         n_pcs = rep(10, 6),
                         base_col = "#7C00FF",
                         base_dir = "analysis/dna_methylation/exploratory/superenhancer_int/")

# Super-enhancers
methyl_se_variances <- variable_features(methyl_se_agg_m)

GenerateExploratoryPlots(data = methyl_se_agg_m,
                         top_features = methyl_se_variances,
                         n_pcs = rep(10, 6),
                         base_col = "#7C00FF",
                         base_dir = "analysis/dna_methylation/exploratory/superenhancer/")

# Enhancers
methyl_enh_variances <- variable_features(methyl_enh_agg_m)

GenerateExploratoryPlots(data = methyl_enh_agg_m,
                         top_features = methyl_enh_variances,
                         n_pcs = rep(10, 6),
                         base_col = "#7C00FF",
                         base_dir = "analysis/dna_methylation/exploratory/enhancer/")

# H3K27ac Peaks
methyl_peak_variances <- variable_features(methyl_peak_agg_m)

GenerateExploratoryPlots(data = methyl_peak_agg_m,
                         top_features = methyl_peak_variances,
                         n_pcs = rep(10, 6),
                         base_col = "#7C00FF",
                         base_dir = "analysis/dna_methylation/exploratory/peak/")

# Gene bodies
methyl_gb_variances <- variable_features(methyl_gb_agg_m)

GenerateExploratoryPlots(data = methyl_gb_agg_m,
                         top_features = methyl_gb_variances,
                         n_pcs = rep(10, 6),
                         base_col = "#7C00FF",
                         base_dir = "analysis/dna_methylation/exploratory/gene_body/")

# Promoters
methyl_tss_variances <- variable_features(methyl_tss_agg_m)

GenerateExploratoryPlots(data = methyl_tss_agg_m,
                         top_features = methyl_tss_variances,
                         n_pcs = rep(10, 6),
                         base_col = "#7C00FF",
                         base_dir = "analysis/dna_methylation/exploratory/promoter/")

# Analyse metabolomics data ----------------------------------------------------------------------------

metab_variances <- variable_features(metab_norm)
metab_ann_col <- data.frame(row.names = colnames(metab_norm),
                            type = factor(metab_metadata$tnbc.subtype))

GenerateExploratoryPlots(data = metab_norm,
                         top_features = metab_variances,
                         n_pcs = rep(10, 6),
                         base_col = "#00CBD4",
                         base_dir = "analysis/metabolomics/exploratory",
                         annotation_col = metab_ann_col)

metab_rep_variances <- variable_features(metab_rep_norm)
metab_rep_ann_col <- data.frame(row.names = colnames(metab_rep_norm), 
                                type = factor(metab_rep_metadata$tnbc.subtype))

GenerateExploratoryPlots(data = metab_rep_norm,
                         top_features = metab_rep_variances,
                         n_pcs = rep(10, 6),
                         base_col = "#00CBD4",
                         base_dir = "analysis/metabolomics/exploratory_replicates",
                         annotation_col = metab_rep_ann_col)

# Analyse HMS data ----------------------------------------------------------------------------

hms_variances <- variable_features(hms_norm)

GenerateExploratoryPlots(data = hms_norm,
                         top_features = hms_variances,
                         n_pcs = rep(10, 6),
                         base_col = "#7A0014",
                         base_dir = "analysis/hms/exploratory/")

hms_rep_variances <- variable_features(hms_rep_norm)
hms_rep_ann_col <- data.frame(row.names = colnames(hms_rep_norm),
                              type = factor(hms_rep_metadata$tnbc.subtype))

GenerateExploratoryPlots(data = hms_rep_norm,
                         top_features = hms_rep_variances,
                         n_pcs = rep(10, 6),
                         base_col = "#7A0014",
                         base_dir = "analysis/hms/exploratory_replicates/",
                         annotation_col = hms_rep_ann_col)

# Analyse BH3 profiling data ----------------------------------------------------------------------------

bh3_ranks <- variable_features(bh3_conc)

GenerateExploratoryPlots(data = bh3_conc,
                         top_features = bh3_ranks,
                         n_pcs = rep(10, 6),
                         base_col = "#3269a8",
                         base_dir = "analysis/bh3_profiling/exploratory/",
                         cap = FALSE)

# Analyse Drug Screen data ----------------------------------------------------------------------------

ds_data <- 1 - ds_auc
ds_ranks <- variable_features(ds_data)

GenerateExploratoryPlots(data = ds_data,
                         top_features = ds_ranks,
                         n_pcs = rep(10, 6),
                         base_col = "#9ea832",
                         base_dir = "analysis/drug_screen/exploratory/",
                         cap = FALSE)