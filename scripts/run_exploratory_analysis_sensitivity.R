library(ggplot2)
library(ggrepel)
library(umap)
library(pheatmap)

setwd("../")

source("scripts/utilities/exploratory_plots.R")
source("scripts/utilities/preprocessing.R")

# Load  data ----------------------------------------------------------------------------

load(file = "R_Data/Heatmap_Metadata.RData")

load(file = "R_Data/Metab_Alt_Data.RData")

# Analyse metabolomics data ----------------------------------------------------------------------------

metab_variances <- variable_features(metab_norm)
metab_ann_col <- data.frame(row.names = colnames(metab_norm),
                            type = factor(metab_metadata$tnbc.subtype))

GenerateExploratoryPlots(data = metab_norm,
                         top_features = metab_variances,
                         n_pcs = rep(10, 6),
                         base_col = "#00CBD4",
                         base_dir = "analysis/metabolomics/sensitivity/exploratory_log10",
                         annotation_col = metab_ann_col)

# Analyse DNA methylation data ----------------------------------------------------------------------------

load(file = "R_Data/Methyl_Data_Overlap_Sensitivity.RData")

methyl_se_variances <- variable_features(methyl_se_agg_m)

GenerateExploratoryPlots(data = methyl_se_agg_m,
                         top_features = methyl_se_variances,
                         n_pcs = rep(10, 6),
                         base_col = "#7C00FF",
                         base_dir = "analysis/dna_methylation/overlap_sensitivity/exploratory/superenhancer/")

load(file = "R_Data/Methyl_Data_CG_Sensitivity.RData")

methyl_se_variances <- variable_features(methyl_se_agg_m)

GenerateExploratoryPlots(data = methyl_se_agg_m,
                         top_features = methyl_se_variances,
                         n_pcs = rep(10, 6),
                         base_col = "#7C00FF",
                         base_dir = "analysis/dna_methylation/cg_sensitivity/exploratory/superenhancer/")