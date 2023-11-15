options(warn = 1)
library(ggsci)
library(survival)
library(ggplot2)
source("../utilities/utils.R")
source("../utilities/primary_data.R")
source("../utilities/preprocessing.R")
source("../utilities/exploratory_plots.R")

setwd("../")

##-----------------------------------------------------------------------------------------------------------------------
#' Load data
load("R_Data/Metabric_Annotated.RData")
load("R_Data/Heatmap_Metadata.RData")
load("R_Data/PRRX1_Targets.RData")
load("R_Data/Gene_Sigs.RData")
load("R_Data/Gene_Sigs_Metadata.RData")

mbc_tnbc_sig_data <- sweep(mbc_tnbc_rna_med, 1, rowMeans(mbc_tnbc_rna_med))

out_path <- "analysis/human_primary_data/metabric/"
alpha <- 0.05

cl_pal <- pal_d3()
cl_cols <- cl_pal(3)

tnbc_types <- c("mes", "bas", "lum")
plot_pct <- 20

cov_ids1 <- c("age_at_diagnosis", "tumor_size", "num_pos_lymph_nodes")
sl <- c(os = "Overall Survival Probability", 
        dss = "Disease-Specific Survival Probability")

##-----------------------------------------------------------------------------------------------------------------------
#' Analysis

analyse_primary_data_correlations(exp_data = mbc_tnbc_sig_data,
                                  exp_metadata = mbc_tnbc_metadata,
                                  imm_scores = imm_n_sig_scores_mbc,
                                  imm_metadata = imm_n_map,
                                  imm_quantiles = imm_n_sig_quantiles_mbc,
                                  patient_metadata = mbc_tnbc_patient_metadata,
                                  data_id = "metabric",
                                  alpha = 0.05,
                                  n_sum_col = 2,
                                  cov_ids1 = cov_ids1)

# Correlate PRRX1 and TNBC type with cellularity
mbc_tnbc_metadata$cellularity <- mbc_tnbc_patient_metadata$cellularity

keep <- mbc_tnbc_metadata$cellularity %in% c("low", "moderate", "high")
cl_exp_data <- mbc_tnbc_sig_data[, keep]
cl_exp_metadata <- droplevels(mbc_tnbc_metadata[keep, ])
cl_imm_n_scores <- imm_n_sig_scores_mbc[keep, ]
ann_cols <- cl_cols
names(ann_cols) <- c("low", "moderate", "high")

analyse_grouping_by_type(type = cl_exp_metadata$cellularity, 
                         groupings = data.frame(tnbc.type = cl_exp_metadata$tnbc.type),
                         path = paste0(out_path, "cellularity/cellularity_and_tnbc_type/"),
                         category_id = "tnbc_type",
                         adjust_p_vals = FALSE,
                         alpha = alpha)

analyse_signature_by_type(type = cl_exp_metadata$cellularity,
                          sig_scores = data.frame(prrx1.expression = as.numeric(cl_exp_data["PRRX1",])),
                          path = paste0(out_path, "cellularity/cellularity_and_prrx1_exp/"),
                          category_id = "cellularity",
                          adjust_p_vals = FALSE,
                          alpha = alpha,
                          ann_col = ann_cols)

analyse_signature_by_type(type = cl_exp_metadata$cellularity,
                          sig_scores = data.frame(hs578.rna.targets = cl_exp_metadata$hs578_rna_targets),
                          path = paste0(out_path, "cellularity/cellularity_and_hs578_rna_targets/"),
                          category_id = "cellularity",
                          adjust_p_vals = FALSE,
                          alpha = alpha,
                          ann_col = ann_cols)

analyse_signature_by_type(type = cl_exp_metadata$cellularity,
                          sig_scores = data.frame(mes.rna.targets = cl_exp_metadata$mes_rna_targets),
                          path = paste0(out_path, "cellularity/cellularity_and_mes_rna_targets/"),
                          category_id = "cellularity",
                          adjust_p_vals = FALSE,
                          alpha = alpha,
                          ann_col = ann_cols)

# Cellularity and immune sigs
analyse_signature_by_type(type = cl_exp_metadata$cellularity, 
                          sig_scores = cl_imm_n_scores, 
                          path = paste0(out_path, "cellularity/cellularity_and_immune_signatures/"), 
                          category_id = "cellularity",
                          alpha = alpha,
                          ann_col = ann_cols)

# Survival analyses
surv_data <- mbc_tnbc_patient_metadata[colnames(cl_exp_data), ]
sp <- "analysis/human_primary_data/metabric/cellularity/cellularity_survival/"

analyse_survival("os_status",  sl["os"],  "os",  "os_months",  "cellularity", surv_data, sp,
                 cov_ids1 = cov_ids1)
analyse_survival("dss_status", sl["dss"], "dss", "dss_months", "cellularity", surv_data, sp,
                 cov_ids1 = cov_ids1)

surv_data <- cbind(mbc_tnbc_metadata, mbc_tnbc_patient_metadata)
sp <- "analysis/human_primary_data/metabric/prrx1_survival/"

for(i in 1:length(tnbc_types)){
  sid <- tnbc_types[i]
  surv_sub <- surv_data[which(surv_data$tnbc.type == tnbc_types[i]),]
  analyse_survival("os_status",  sl["os"],  "os",  "os_months",  "mes_rna_targets",  surv_sub, sp, 
                   cov_ids1 = cov_ids1, samples_id = sid)
  analyse_survival("dss_status", sl["dss"], "dss", "dss_months",  "mes_rna_targets", surv_sub, sp, 
                   cov_ids1 = cov_ids1, samples_id = sid)
  
  analyse_survival("os_status",  sl["os"],  "os",  "os_months",  "hs578_rna_targets", surv_sub, sp, 
                   cov_ids1 = cov_ids1, samples_id = sid)
  analyse_survival("dss_status", sl["dss"], "dss", "dss_months", "hs578_rna_targets", surv_sub, sp, 
                   cov_ids1 = cov_ids1, samples_id = sid)
  
  analyse_survival("os_status",  sl["os"],  "os",  "os_months",  "prrx1", surv_sub, sp, 
                   cov_ids1 = cov_ids1, samples_id = sid)
  analyse_survival("dss_status", sl["dss"], "dss", "dss_months", "prrx1", surv_sub, sp, 
                   cov_ids1 = cov_ids1, samples_id = sid)
}

##-----------------------------------------------------------------------------------------------------------------------
#' Exploratory analysis

plot_path <- create_path(paste0(out_path, "exploratory/"))

var_feats <- variable_features(mbc_tnbc_rna_med)
ntop <- round((plot_pct * length(var_feats)) / 100)
ann_col <- data.frame(row.names = rownames(mbc_tnbc_metadata), type = mbc_tnbc_metadata$tnbc.type)

pdf(file.path(plot_path, "tnbc_heatmap.pdf"))
res <- PlotHeatmap(mbc_tnbc_rna_med[var_feats[1:ntop], ], 
                   annotation_col = ann_col,
                   annotation_colours = list(type = line_colours$type2),
                   label_samples = FALSE)
dev.off()

mbc_tnbc_genes <- unique(unlist(uniform_subset_gsl(gsl = tnbc_gsl, genes = rownames(mbc_tnbc_rna_med))))

pdf(file.path(plot_path, "tnbc_heatmap_unif_type_genes.pdf"))
res <- PlotHeatmap(mbc_tnbc_rna_med[mbc_tnbc_genes, ], 
                   annotation_col = ann_col,
                   annotation_colours = list(type = line_colours$type2),
                   label_samples = FALSE)
dev.off()

mbc_tnbc_genes <- unique(unlist(subset_gsl(gsl = tnbc_gsl, genes = rownames(mbc_tnbc_rna_med))))

pdf(file.path(plot_path, "tnbc_heatmap_type_genes.pdf"))
res <- PlotHeatmap(mbc_tnbc_rna_med[mbc_tnbc_genes, ], 
                   annotation_col = ann_col,
                   annotation_colours = list(type = line_colours$type2))
dev.off()

##-----------------------------------------------------------------------------------------------------------------------
#' Summary stats for writeup

median(mbc_tnbc_patient_metadata$os_months)
median(mbc_tnbc_patient_metadata$dss_months)
