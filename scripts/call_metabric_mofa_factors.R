#####
# assign factor scores to Metabric TNBC samples using RNA data
#####

setwd("../")

################################################################################
# load library
################################################################################

library(MOFA2)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(mclust)
source("../utilities/utils.R")
source("../utilities/plots.R")
source("../utilities/mofa2_factor_calls.R")

################################################################################
# load data
################################################################################

load("R_Data/Metabric_Annotated.RData")
load("R_Data/Heatmap_Metadata.RData")
mofa_object <- readRDS("R_Data/MOFA2_TopN_Heatmap_Feature_Filtering_Seed20200225.hdf5")
mofa_weights <- get_weights(mofa_object)
out_path <- "analysis/human_primary_data/metabric/"
factor_assocs <- c(factor_2 = "lum", factor_3 = "mes", factor_6 = "bas")

################################################################################
# transform data
################################################################################

# create rna dataset list
rna_data_list <- list(rna = mbc_tnbc_rna_med)

rna_data_list <- format_model_data_complete_samples(rna_data_list)

#' prepare data to call factors:
#' combine expression and weight data lists into a single object, center
#' expression features, subset and reorder both expression and weight features
#' to match
rna_weight_list <- mofa_weights[names(rna_data_list)]
rna_mofa_call_data <- prepare_mofa_calls(rna_data_list, rna_weight_list)

################################################################################
# local functions
################################################################################

plot_variance_explained_mofa_calls <- function(gg_ve, pal = ds_colours$dataset){
  p <- ggplot(gg_ve, aes(x = dataset, y = mbc_tnbc_r2, fill = dataset)) +
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
rna_mofa_scores <- call_mofa_factor_scores(rna_mofa_call_data)

# quantize factor scores
rna_mofa_quantiles <- quantize(rna_mofa_scores)

# Analyse variance explained in each dataset
rna_var_exp <- get_variance_explained_mofa_calls(rna_mofa_scores, rna_mofa_call_data)

gg_r_ve <- data.frame(dataset = names(rna_var_exp$r2_total),
                      mbc_tnbc_r2 = rna_var_exp$r2_total)
pdf(file.path(out_dir, "mbc_tnbc_rna_factors_r2.pdf"))
plot_variance_explained_mofa_calls(gg_r_ve)
dev.off()

pdf(file.path(out_dir, "mbc_tnbc_rna_factors.pdf"))
par(mfrow = c(2, 4))
tmp <- lapply(colnames(rna_mofa_scores), function(i){
  x <- rna_mofa_scores[,i]
  hist(x, breaks = 50, main = i, xlab = "Factor Score")
})
dev.off()

# rna
mbc_mofa_factors <- as.data.frame(rna_mofa_scores)
mbc_mofa_factor_quantiles <- rna_mofa_quantiles

################################################################################
# Annotate factor data
################################################################################

clust_path <- create_path(file.path(out_dir, "clustering"))

# create metadata
metadata_cols <- c("tnbc.type", "prrx1", "hs578_rna_targets", "mes_rna_targets")
all(rownames(mbc_mofa_factors) %in% rownames(mbc_tnbc_metadata))
mbc_mofa_metadata <- mbc_tnbc_metadata[rownames(mbc_mofa_factors), metadata_cols]

all(rownames(mbc_mofa_metadata) %in% rownames(imm_n_sig_scores_mbc))
imm_n_sig_scores_mbc_mofa <- imm_n_sig_scores_mbc[rownames(mbc_mofa_metadata),]

all(rownames(mbc_mofa_metadata) %in% rownames(mbc_tnbc_patient_metadata))
mbc_mofa_patient_metadata <- mbc_tnbc_patient_metadata[rownames(mbc_mofa_metadata), ]

# assign factor 3 clusters
fv <- mbc_mofa_factors[, c("Factor3")]
mclust_res <- Mclust(fv)
pdf(file.path(clust_path, "mclust_refined_type.pdf"))
plot_mclust_results(mclust_res)
dev.off()

stopifnot(length(unique(mclust_res$classification)) == 2)
clust_means <- tapply(fv, mclust_res$classification, mean)
stopifnot(clust_means["2"] > clust_means["1"])
stopifnot(all(mbc_mofa_metadata$tnbc.type[which(mclust_res$classification == 2)] == "mes"))
mbc_mofa_metadata$refined.type <- sapply(1:nrow(mbc_mofa_metadata), function(x){
  if(mclust_res$classification[x] == 2){
    "mes_high"
  } else if (mbc_mofa_metadata$tnbc.type[x] == "mes"){
    "mes_low"
  } else {
    as.character(mbc_mofa_metadata$tnbc.type[x])
  }
})
mbc_mofa_metadata$refined.type <- factor(mbc_mofa_metadata$refined.type,
                                         levels = c("bas", "lum", "mes_low", "mes_high"))

# simple threshold
f3_threshold <- max(mbc_mofa_factors$Factor3[which(mbc_mofa_metadata$tnbc.type != "mes")])
mbc_mofa_metadata$factor.3.group <- sapply(1:nrow(mbc_mofa_metadata), function(x){
  if(mbc_mofa_factors$Factor3[x] > f3_threshold){
    "mes_high"
  } else if (mbc_mofa_metadata$tnbc.type[x] == "mes"){
    "mes_low"
  } else {
    as.character(mbc_mofa_metadata$tnbc.type[x])
  }
})
mbc_mofa_metadata$factor.3.group <- factor(mbc_mofa_metadata$factor.3.group,
                                           levels = c("bas", "lum", "mes_low", "mes_high"))

# assign mofa types
fv <- mbc_mofa_factors[, c("Factor2", "Factor3", "Factor6")]
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
mbc_mofa_metadata$mofa.type <- factor(easy_find_replace(as.character(mclust_res$classification), 
                                                        c("1", "2", "3"),
                                                        clust_names),
                                      levels = c("mofa_bas", "mofa_lum", "mofa_mes"))

################################################################################
# save results
################################################################################

save(mbc_mofa_factors, 
     mbc_mofa_factor_quantiles, 
     mbc_mofa_metadata,
     mbc_mofa_patient_metadata,
     imm_n_sig_scores_mbc_mofa,
     file = "R_Data/Metabric_MOFA_Data.RData")
save(rna_mofa_scores, rna_var_exp,file = "R_Data/Metabric_MOFA_Intermediate_Data.RData")
