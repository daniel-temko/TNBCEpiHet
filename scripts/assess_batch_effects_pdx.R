library(ggsci)

setwd("../")

source("scripts/utilities/exploratory_plots.R")
source("scripts/utilities/preprocessing.R")

# Load  data ----------------------------------------------------------------------------

load("R_Data/PDX_Methyl_Data.RData")

batch_ids <- c("Batch_1", "Batch_2")

pdx_methyl_metadata$batch <- batch_ids[match(pdx_methyl_metadata$Date.Arrayed, c("20170127", "20170620"))]

batch_pal <- pal_d3()
batch_cols <- batch_pal(2)

# Analyse DNA methylation data ----------------------------------------------------------------------------

pdx_methyl_se_agg_variances <- variable_features(pdx_methyl_se_agg_m)
pdx_methyl_ann_col <- data.frame(row.names = colnames(pdx_methyl_se_agg_m),
                                 type = factor(pdx_methyl_metadata$batch))
pdx_methyl_ann_colours <- list(type = c(Batch_1 = batch_cols[1], Batch_2 = batch_cols[2]))

GenerateExploratoryPlots(data = pdx_methyl_se_agg_m,
                         top_features = pdx_methyl_se_agg_variances,
                         n_pcs = rep(10, 6),
                         base_col = "#7C00FF",
                         base_dir = "analysis/pdx_dna_methylation/qc/exploratory_by_batch/superenhancer",
                         annotation_col = pdx_methyl_ann_col,
                         annotation_colours = pdx_methyl_ann_colours)
