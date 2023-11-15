#####
# assign factor scores to TCGA TNBC samples using RNA and DNA methylation data
#####

setwd("../")

################################################################################
# load library
################################################################################

library(MOFA2)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(mclust)
source("../utilities/utils.R")
source("../utilities/plots.R")
source("../utilities/mofa2_factor_calls.R")

################################################################################
# load data
################################################################################

load("R_Data/TCGA_RNA_Annotated.RData")
load("R_Data/TCGA_Methyl_Data.RData")
load("R_Data/Heatmap_Metadata.RData")
mofa_object <- readRDS("R_Data/MOFA2_TopN_Heatmap_Feature_Filtering_Seed20200225.hdf5")
mofa_weights <- get_weights(mofa_object)
out_path <- "analysis/human_primary_data/tcga"
factor_assocs <- c(factor_2 = "lum", factor_3 = "mes", factor_6 = "bas")

################################################################################
# transform data
################################################################################

# Update methylation feature names
rownames(tcga_tnbc_methyl_cl_se_int_agg_m) <- paste0(rownames(tcga_tnbc_methyl_cl_se_int_agg_m), "_methyl")
rownames(tcga_tnbc_methyl_gb_agg_m) <- paste0(rownames(tcga_tnbc_methyl_gb_agg_m), "_gb_methyl")
rownames(tcga_tnbc_methyl_tss_agg_m) <- paste0(rownames(tcga_tnbc_methyl_tss_agg_m), "_tss_methyl")

# create methyl_rna dataset list
methyl_rna_data_list <- list(methyl.se = tcga_tnbc_methyl_cl_se_int_agg_m,
                             methyl.gb = tcga_tnbc_methyl_gb_agg_m,
                             methyl.tss = tcga_tnbc_methyl_tss_agg_m,
                             rna = tcga_tnbc_rna_fpkm_norm_uf)

#' Remove samples missing for any dataset and convert to matrix with matching
#' sample orders across datasets
methyl_rna_data_list <- format_model_data_complete_samples(methyl_rna_data_list)

#' prepare data to call factors:
#' combine expression and weight data lists into a single object, center
#' expression features, subset and reorder both expression and weight features
#' to match
methyl_rna_weight_list <- mofa_weights[names(methyl_rna_data_list)]
methyl_rna_mofa_call_data <- prepare_mofa_calls(methyl_rna_data_list, methyl_rna_weight_list)

################################################################################
# local functions
################################################################################

plot_variance_explained_mofa_calls <- function(gg_ve, pal = ds_colours$dataset){
  p <- ggplot(gg_ve, aes(x = dataset, y = tcga_tnbc_r2, fill = dataset)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = pal) +
    ylim(c(0, 0.5)) + 
    my_theme()
  return(p)
}

quantize <- function(mofa_scores){
  mofa_quantiles <- as.data.frame(lapply(1:ncol(mofa_scores), function(x) {
    get_quantiles(mofa_scores[,x], 2)
  }))
  rownames(mofa_quantiles) <- rownames(mofa_scores)
  colnames(mofa_quantiles) <- paste0(colnames(mofa_scores), ".quantile")
  return(mofa_quantiles)
}

################################################################################
# call factors
################################################################################

# output directory for plots
out_dir <- file.path(out_path, "mofa_factors")
if(!dir.exists(out_dir)) dir.create(out_dir)

# call factors
methyl_rna_mofa_scores <- call_mofa_factor_scores(methyl_rna_mofa_call_data)

# quantize factor scores
methyl_rna_mofa_quantiles <- quantize(methyl_rna_mofa_scores)

# Analyze variance explained in each dataset
methyl_rna_var_exp <- get_variance_explained_mofa_calls(methyl_rna_mofa_scores, methyl_rna_mofa_call_data)

gg_mr_ve <- data.frame(dataset = names(methyl_rna_var_exp$r2_total),
                       tcga_tnbc_r2 = methyl_rna_var_exp$r2_total)
pdf(file.path(out_dir, "tcga_tnbc_methyl_rna_factors_r2.pdf"))
plot_variance_explained_mofa_calls(gg_mr_ve)
dev.off()

pdf(file.path(out_dir, "tcga_tnbc_methyl_rna_factors.pdf"))
par(mfrow = c(2, 4))
tmp <- lapply(colnames(methyl_rna_mofa_scores), function(i){
  x <- methyl_rna_mofa_scores[,i]
  hist(x, breaks = 50, main = i, xlab = "Factor Score")
})
dev.off()

# methyl_rna
tcga_mofa_factors <- as.data.frame(methyl_rna_mofa_scores)
tcga_mofa_factor_quantiles <- methyl_rna_mofa_quantiles

################################################################################
# Annotate factor data
################################################################################

clust_path <- create_path(file.path(out_dir, "clustering"))

metadata_cols <- c("tnbc.type", "prrx1", "hs578_rna_targets", "mes_rna_targets")
all(rownames(tcga_mofa_factors) %in% rownames(tcga_tnbc_rna_metadata))
tcga_mofa_metadata <- tcga_tnbc_rna_metadata[rownames(tcga_mofa_factors), metadata_cols]

all(rownames(tcga_mofa_metadata) %in% rownames(imm_n_sig_scores_tcga))
imm_n_sig_scores_tcga_mofa <- imm_n_sig_scores_tcga[rownames(tcga_mofa_metadata), ]

all(rownames(tcga_mofa_metadata) %in% rownames(tcga_tnbc_rna_patient_metadata))
tcga_mofa_patient_metadata <- tcga_tnbc_rna_patient_metadata[rownames(tcga_mofa_metadata), ]

# assign factor 3 clusters
fv <- tcga_mofa_factors[, "Factor3"]
mclust_res <- Mclust(fv)
pdf(file.path(clust_path, "mclust_refined_type.pdf"))
plot_mclust_results(mclust_res)
dev.off()

stopifnot(length(unique(mclust_res$classification)) == 2)
clust_means <- tapply(fv, mclust_res$classification, mean)
stopifnot(clust_means["2"] > clust_means["1"])
stopifnot(all(tcga_mofa_metadata$tnbc.type[which(mclust_res$classification == 2)] == "mes"))
tcga_mofa_metadata$refined.type <- sapply(1:nrow(tcga_mofa_metadata), function(x){
  if(mclust_res$classification[x] == 2){
    "mes_high"
  } else if (tcga_mofa_metadata$tnbc.type[x] == "mes"){
    "mes_low"
  } else {
    as.character(tcga_mofa_metadata$tnbc.type[x])
  }
})
tcga_mofa_metadata$refined.type <- factor(tcga_mofa_metadata$refined.type,
                                          levels = c("bas", "lum", "mes_low", "mes_high"))

# simple threshold
f3_threshold <- max(tcga_mofa_factors$Factor3[which(tcga_mofa_metadata$tnbc.type != "mes")])
tcga_mofa_metadata$factor.3.group <- sapply(1:nrow(tcga_mofa_metadata), function(x){
  if(tcga_mofa_factors$Factor3[x] > f3_threshold){
    "mes_high"
  } else if (tcga_mofa_metadata$tnbc.type[x] == "mes") {
    "mes_low"
  } else {
    as.character(tcga_mofa_metadata$tnbc.type[x])
  }
})
tcga_mofa_metadata$factor.3.group <- factor(tcga_mofa_metadata$factor.3.group,
                                            levels = c("bas", "lum", "mes_low", "mes_high"))


# assign mofa types
fv <- tcga_mofa_factors[, c("Factor2", "Factor3", "Factor6")]
mclust_res <- Mclust(fv)
pdf(file.path(clust_path, "mclust_mofa_type.pdf"))
plot_mclust_results(mclust_res)
dev.off()

stopifnot(length(unique(mclust_res$classification)) == 3)
clust_means <- apply(fv, 2, function(x){
  tapply(x, mclust_res$classification, mean)
})
old_means <- apply(tcga_mofa_factors[, c("Factor2", "Factor3", "Factor6")], 2, function(x){
  tapply(x, tcga_mofa_metadata$tnbc.type, mean)
})
dist(rbind(clust_means, old_means))
factor_clusts <- apply(clust_means, 2, which.max)
stopifnot(length(unique(factor_clusts)) == 3)
clust_names <- paste0("mofa_", factor_assocs[factor_clusts])
tcga_mofa_metadata$mofa.type <- factor(easy_find_replace(as.character(mclust_res$classification),
                                                         c("1", "2", "3"),
                                                         clust_names),
                                       levels = c("mofa_bas", "mofa_lum", "mofa_mes"))


################################################################################
# save results
################################################################################

save(tcga_mofa_factors, 
     tcga_mofa_factor_quantiles,
     tcga_mofa_metadata,
     tcga_mofa_patient_metadata,
     imm_n_sig_scores_tcga_mofa,
     file = "R_Data/TCGA_MOFA_Data.RData")
save(methyl_rna_mofa_scores,
     methyl_rna_var_exp,
     file = "R_Data/TCGA_MOFA_Intermediate_Data.RData")
