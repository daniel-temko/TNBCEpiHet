library(ggplot2)
library(ggrepel)
library(umap)
library(ggsci)
library(pheatmap)

setwd("../")

source("scripts/utilities/exploratory_plots.R")
source("scripts/utilities/preprocessing.R")

# Load  data ----------------------------------------------------------------------------

load(file = "R_Data/RNA_Data.RData")
load(file = "R_Data/Methyl_Data.RData")
load(file = "R_Data/Metab_Data.RData")

batch_ids <- c("Batch_1", "Batch_2")

rna_metadata$batch <- batch_ids[match(rna_metadata$Date.Sequenced, c("20160104", "20160918"))]
methyl_metadata$batch <- batch_ids[match(methyl_metadata$Date.Arrayed, c("20161028", "20170127"))]
metab_rep_metadata$batch2 <- paste0("Batch_", metab_rep_metadata$batch)

batch_pal <- pal_d3()
batch_cols <- batch_pal(2)

# Analyse RNA-seq data ----------------------------------------------------------------------------

rna_ranks <- variable_features(rna_rpkm_norm)
rna_ann_col <- data.frame(row.names = colnames(rna_rpkm_norm),
                          type = factor(rna_metadata$batch))
rna_ann_colours <- list(type = c(Batch_1 = batch_cols[1], Batch_2 = batch_cols[2]))

GenerateExploratoryPlots(data = rna_rpkm_norm,
                         top_features = rna_ranks,
                         n_pcs = rep(10, 6),
                         base_col = "#FF00FF",
                         base_dir = "analysis/rna_seq/qc/exploratory_by_batch",
                         annotation_col = rna_ann_col,
                         annotation_colours = rna_ann_colours)

# Analyse DNA methylation data ----------------------------------------------------------------------------

methyl_se_agg_variances <- variable_features(methyl_se_agg_m)
methyl_ann_col <- data.frame(row.names = colnames(methyl_se_agg_m),
                             type = factor(methyl_metadata$batch))
methyl_ann_colours <- list(type = c(Batch_1 = batch_cols[1], Batch_2 = batch_cols[2]))

GenerateExploratoryPlots(data = methyl_se_agg_m,
                         top_features = methyl_se_agg_variances,
                         n_pcs = rep(10, 6),
                         base_col = "#7C00FF",
                         base_dir = "analysis/dna_methylation/qc/exploratory_by_batch/superenhancer",
                         annotation_col = methyl_ann_col,
                         annotation_colours = methyl_ann_colours)

# Analyse metabolomics data ----------------------------------------------------------------------------

metab_rep_variances <- variable_features(metab_rep_norm)
metab_rep_ann_col <- data.frame(row.names = colnames(metab_rep_norm),
                                type = factor(metab_rep_metadata$batch2))
metab_rep_ann_colours <- list(type = c(Batch_1 = batch_cols[1], Batch_2 = batch_cols[2]))

GenerateExploratoryPlots(data = metab_rep_norm,
                         top_features = metab_rep_variances,
                         n_pcs = rep(10, 6),
                         base_col = "#00CBD4",
                         base_dir = "analysis/metabolomics/qc/exploratory_by_batch",
                         annotation_col = metab_rep_ann_col,
                         annotation_colours = metab_rep_ann_colours)
