options(warn = 2)
source("../utilities/utils.R")
source("../utilities/primary_data.R")
source("../utilities/preprocessing.R")
source("../utilities/exploratory_plots.R")

setwd("../")

##-----------------------------------------------------------------------------------------------------------------------
#' Load data
load("R_Data/TCGA_RNA_Annotated.RData")
load("R_Data/TCGA_TNBC_CibersortX.RData")
load("R_Data/Heatmap_Metadata.RData")
load("R_Data/PRRX1_Targets.RData")
load("R_Data/Gene_Sigs.RData")
load("R_Data/Gene_Sigs_Metadata.RData")

out_path <- "analysis/human_primary_data/tcga/"

tcga_tnbc_sig_data <- sweep(tcga_tnbc_rna_fpkm_norm_uf, 1, rowMeans(tcga_tnbc_rna_fpkm_norm_uf))

plot_pct <- 20

cov_ids1 <- c("age", "path_stage")
#cov_ids2 <- c("age", "path_stage", "radiation")

##-----------------------------------------------------------------------------------------------------------------------
#' Analysis

analyse_primary_data_correlations(exp_data = tcga_tnbc_sig_data,
                                  exp_metadata = tcga_tnbc_rna_metadata,
                                  imm_scores = imm_n_sig_scores_tcga,
                                  imm_metadata = imm_n_map,
                                  imm_quantiles = imm_n_sig_quantiles_tcga,
                                  patient_metadata = tcga_tnbc_rna_patient_metadata,
                                  data_id = "tcga",
                                  alpha = 0.05,
                                  n_sum_col = 3,
                                  cov_ids1 = cov_ids1,
                                  test_pfs = TRUE)

# Correlate TNBC type with CIBERSORTx absolute signature score
analyse_signature_by_type(type = tcga_tnbc_rna_metadata$tnbc.type,
                          sig_scores = data.frame(Total.Absolute = tcga_tnbc_cibersortx$Absolute.score..sig.score.),
                          path = paste0(out_path, "cibersortx/absolute_sig_and_tnbc_type/"),
                          category_id = "tnbc_type",
                          adjust_p_vals = FALSE,
                          alpha = 0.05,
                          ann_col = line_colours$type2)

##-----------------------------------------------------------------------------------------------------------------------
#' Exploratory analysis

plot_path <- create_path(paste0(out_path, "exploratory/"))

# filter genes
keep <- apply(tcga_tnbc_rna_fpkm_norm_uf, 1, var) > 0
tcga_exp_data <- tcga_tnbc_rna_fpkm_norm_uf[keep, ]

keep <- apply(tcga_exp_data, 1, function(x) length(which(x > 1))) > 1
tcga_exp_data <- tcga_exp_data[keep, ]

var_feats <- variable_features(tcga_exp_data)
ntop <- ((plot_pct * length(var_feats)) / 100)

ann_col <- data.frame(row.names = rownames(tcga_tnbc_rna_metadata), type = tcga_tnbc_rna_metadata$tnbc.type)

pdf(file.path(plot_path, "tnbc_heatmap.pdf"))
res <- PlotHeatmap(tcga_exp_data[var_feats[1:ntop], ],
                   annotation_col = ann_col,
                   annotation_colours = list(type = line_colours$type2),
                   label_samples = FALSE)
dev.off()

tcga_tnbc_genes <- unique(unlist(uniform_subset_gsl(gsl = tnbc_gsl, rownames(tcga_tnbc_rna_fpkm_norm_uf))))

pdf(file.path(plot_path, "tnbc_heatmap_unif_type_genes.pdf"))
res <- PlotHeatmap(tcga_tnbc_rna_fpkm_norm_uf[tcga_tnbc_genes, ],
                   annotation_col = ann_col,
                   annotation_colours = list(type = line_colours$type2),
                   label_samples = FALSE)
dev.off()

tcga_tnbc_genes <- unique(unlist(subset_gsl(gsl = tnbc_gsl, rownames(tcga_tnbc_rna_fpkm_norm_uf))))

pdf(file.path(plot_path, "tnbc_heatmap_type_genes.pdf"))
res <- PlotHeatmap(tcga_tnbc_rna_fpkm_norm_uf[tcga_tnbc_genes, ],
                   annotation_col = ann_col,
                   annotation_colours = list(type = line_colours$type2),
                   label_samples = FALSE)
dev.off()

##-----------------------------------------------------------------------------------------------------------------------
#' Summary stats for writeup

median(tcga_tnbc_rna_patient_metadata$os_months, na.rm = TRUE)
median(tcga_tnbc_rna_patient_metadata$dss_months, na.rm = TRUE)
median(tcga_tnbc_rna_patient_metadata$pfs_months, na.rm = TRUE)
