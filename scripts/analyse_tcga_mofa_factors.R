#####
# Model TCGA survival dependence on MOFA factors
#####

setwd("../")

################################################################################
# load library
################################################################################

source("../utilities/primary_data.R")
source("../utilities/utils.R")
library(cowplot)
library(ggpubr)
library(reshape)

################################################################################
# load data
################################################################################

load(file = "R_Data/Heatmap_Metadata.RData")
load(file = "R_Data/TCGA_MOFA_Data.RData")
load(file = "R_Data/Gene_Sigs_Metadata.RData")
out_path <- "analysis/human_primary_data/tcga/mofa_factors"
cov_ids1 <- c("age", "path_stage")
sl <- c(os = "Overall Survival Probability", 
        dss = "Disease-Specific Survival Probability",
        pfs = "Progression-Free Survival Probability")
refined_type_colors <- c(line_colours$type2[c("bas", "lum")], 
                         "mes_low" = "lightgreen", "mes_high" = "darkgreen")
mofa_type_colors <- line_colours$type2
names(mofa_type_colors) <- paste0("mofa_", names(mofa_type_colors))

alpha <- 0.05

################################################################################
# Local functions
################################################################################

plot_factor_raster <- function(factor_scores, factor_metadata, factor_ids, col_by = "tnbc.type",
                               pal = line_colours$type2){
  stopifnot(is.data.frame(factor_scores) & is.data.frame(factor_metadata))
  plot_list <- list()
  counter <- 1
  for(i in factor_ids){
    for(j in factor_ids){
      xlims <- c(min(factor_scores[[i]]), max(factor_scores[[i]]))
      ylims <- c(min(factor_scores[[j]]), max(factor_scores[[j]]))
      gg <- data.frame(x = factor_scores[[i]],
                       y = factor_scores[[j]],
                       type = factor_metadata[[col_by]])
      p <- ggplot(gg, aes(x = x, y = y, color = type)) + 
        geom_point() + scale_color_manual(values = pal) +
        xlab(i) + ylab(j) + 
        scale_x_continuous(breaks = seq(ceiling(xlims[1] / 0.5) * 0.5, floor(xlims[2] / 0.5) * 0.5, 0.5)) +
        scale_y_continuous(breaks = seq(ceiling(ylims[1] / 0.5) * 0.5, floor(ylims[2] / 0.5) * 0.5, 0.5)) +
        my_theme() 
      if(i != j){
        p <- p + theme(legend.position = "none")
      } else if(counter == 1) {
        p <- as_ggplot(get_legend(p))
      } else{
        p <- NULL
      }
      plot_list[[counter]] <- p
      counter <- counter + 1
    }
  }
  plot_grid(plotlist = plot_list, ncol = length(factor_ids), byrow = FALSE)
}

analyse_signature_by_focused_type <- function(type, type1, type2, sig_scores, path, category_id, adjust_p_vals = TRUE,
                                              alpha = 0.05, ann_col = refined_type_colors, text_size = 30,
                                              lab_size = 8, y_lab_vals = NULL, summary_list = list()){
  plot_path <- create_path(paste0(path, "plots/"))
  list_path <- create_path(paste0(path, "lists/"))
  
  for(i in 1:ncol(sig_scores)){
    sig_id <- colnames(sig_scores)[i]
    message(sig_id)
    
    if(!is.null(y_lab_vals)){
      y_lab <- y_lab_vals[i]
    } else {
      y_lab <- sig_id
    }
    
    pdf(paste0(plot_path, sig_id, "_by_", category_id, ".pdf"), width = 11, height = 9)
    if(min(table(type)) >= 5){
      p <- plot_vioplot(type = type, y = sig_scores[,i], ylab = y_lab,
                        ann_col = ann_col, text_size = text_size, lab_size = lab_size)  
    } else {
      p <- plot_boxplot(type = type, y = sig_scores[,i], ylab = y_lab,
                        ann_col = ann_col, text_size = text_size, lab_size = lab_size)
    }
    print(p)
    dev.off()
    
    # create summary boxplots
    pdf(paste0(plot_path, sig_id, "_by_", category_id, "_compact.pdf"), width = 6, height = 3)
    pp <- data.frame(value = sig_scores[,i], type = type)
    p <- ggplot(pp, aes(x = type, y = value, fill = type)) + 
      geom_boxplot(outlier.shape=NA, show.legend = FALSE) +
      scale_fill_manual(values = ann_col) + 
      geom_jitter(color = "black", size = 2, alpha = 0.9, show.legend = FALSE) +
      my_theme() + theme(text = element_text(size = 24), axis.text.x = element_text(size = 24)) +
      xlab("") + ylab(y_lab)
    print(p)
    dev.off()
  }
  
  # test differences between type 1 and 2 only
  stopifnot(all(c(type1, type2) %in% type))
  keep <- type %in% c(type1, type2)
  sub_type <- type[keep]
  sub_sig <- sig_scores[keep, ]
  
  type1_med <- apply(sub_sig, 2, function(x){
    median(x[which(sub_type == type1)])
  })
  type2_med <- apply(sub_sig, 2, function(x){
    median(x[which(sub_type == type2)])
  })
  p_vals <- apply(sub_sig, 2, function(x){
    test_df <- data.frame(val = x, g = sub_type)
    wilcox.test(x ~ g, test_df)$p.value
  })
  sig_tests <- data.frame(diff.in.medians = type2_med - type1_med, p.val = p_vals)
  
  if(adjust_p_vals){
    sig_tests$p.adj <- p.adjust(p_vals, method = "holm")
  }
  
  sig_tests <- sig_tests[order(sig_tests$p.val), ]
  
  write.table(sig_tests, file = paste0(list_path, category_id, "_wilcox_test_p_values.csv"), 
              sep = ",", col.names = NA)
  
  # create summary heatmap
  for(plot_id in names(summary_list)){
    aa <- aggregate(sig_scores, by = list(type), median)
    at <- t(aa[,-1])
    colnames(at) <- aa[,1]
    at <- at[summary_list[[plot_id]], , drop = FALSE]
    at <- at[order(rownames(at)), , drop = FALSE]
    pdf(paste0(path, plot_id, "_summary.pdf"), width = 6, height = 3)
    p <- pheatmap(at, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 24)
    print(p)
    dev.off() 
  }
}

################################################################################
# Analyze factor call assocations
################################################################################

# analyze associations with TNBC type assignments for type-linked factors
tnbc_type_path <- create_path(file.path(out_path, "mofa_factors_and_tnbc_type"))
refined_type_path <- create_path(file.path(out_path, "mofa_factors_and_refined_types"))
factor_3_group_path <- create_path(file.path(out_path, "mofa_factors_and_factor_3_group"))
mofa_type_path <- create_path(file.path(out_path, "mofa_factors_and_mofa_type"))

# factor 2 and 3 plot by TNBC type
gg <- cbind(tcga_mofa_factors, tnbc.type = tcga_mofa_metadata$tnbc.type)
pdf(file.path(tnbc_type_path, "factor_2_and_3_by_tnbc_type.pdf"))
ggplot(gg, aes(x = Factor2, y = Factor3, color = tnbc.type)) +
  geom_point() +
  scale_color_manual(values = line_colours$type2) +
  my_theme()
dev.off()

# factor 2, 3, and 6 raster by TNBC type
pdf(file.path(tnbc_type_path, "tnbc_type_factors_by_tnbc_type.pdf"), 9, 9)
plot_factor_raster(tcga_mofa_factors, tcga_mofa_metadata, c("Factor2", "Factor3", "Factor6"))
dev.off()

# factor 2 and 3 plot by refined types
gg <- cbind(tcga_mofa_factors, refined.type = tcga_mofa_metadata$refined.type)
pdf(file.path(refined_type_path, "factor_2_and_3_by_refined_type.pdf"))
ggplot(gg, aes(x = Factor2, y = Factor3, color = refined.type)) +
  geom_point() +
  scale_color_manual(values = refined_type_colors) +
  my_theme()
dev.off()

# factor 2, 3, and 6 raster by refined type
pdf(file.path(refined_type_path, "tnbc_type_factors_by_type.pdf"), 9, 9)
plot_factor_raster(tcga_mofa_factors, tcga_mofa_metadata, c("Factor2", "Factor3", "Factor6"),
                   col_by = "refined.type", pal = refined_type_colors)
dev.off()

# factor 2 and 3 plot by factor 3 group
gg <- cbind(tcga_mofa_factors, factor.3.group = tcga_mofa_metadata$factor.3.group)
pdf(file.path(factor_3_group_path, "factor_2_and_3_by_factor_3_group.pdf"))
ggplot(gg, aes(x = Factor2, y = Factor3, color = factor.3.group)) +
  geom_point() +
  scale_color_manual(values = refined_type_colors) +
  my_theme()
dev.off()

# factor 2, 3, and 6 raster by factor 3 group
pdf(file.path(factor_3_group_path, "tnbc_type_factors_by_factor_3_group.pdf"), 9, 9)
plot_factor_raster(tcga_mofa_factors, tcga_mofa_metadata, c("Factor2", "Factor3", "Factor6"),
                   col_by = "factor.3.group", pal = refined_type_colors)
dev.off()

# factor 2 and 3 plot by MOFA type
gg <- cbind(tcga_mofa_factors, mofa.type = tcga_mofa_metadata$mofa.type)
pdf(file.path(mofa_type_path, "factor_2_and_3_by_mofa_type.pdf"))
ggplot(gg, aes(x = Factor2, y = Factor3, color = mofa.type)) +
  geom_point() +
  scale_color_manual(values = mofa_type_colors) +
  my_theme()
dev.off()

# factor 2, 3, and 6 raster by MOFA type
pdf(file.path(mofa_type_path, "tnbc_type_factors_by_mofa_type.pdf"), 9, 9)
plot_factor_raster(tcga_mofa_factors, tcga_mofa_metadata, c("Factor2", "Factor3", "Factor6"),
                   col_by = "mofa.type", pal = mofa_type_colors)
dev.off()

# correspondence between factor 3 group and refined type and covariates

# associations with PRRX1 expresion and target signatures
comp_ids <- c("prrx1", "mes_rna_targets", "hs578_rna_targets")
summary_list <- list(prrx1 = comp_ids[1], targets = comp_ids[2:3])
y_lab_vals <- c("PRRX1 Expr.", "MRT Sig.", "HsRT Sig.")
analyse_signature_by_focused_type(type = tcga_mofa_metadata$refined.type,
                                  type1 = "mes_low",
                                  type2 = "mes_high",
                                  sig_scores = tcga_mofa_metadata[, comp_ids],
                                  path = file.path(out_path, "prrx1_and_refined_type/"),
                                  category_id = "refined_type",
                                  adjust_p_vals = FALSE,
                                  alpha = alpha,
                                  y_lab_vals = y_lab_vals)

analyse_signature_by_focused_type(type = tcga_mofa_metadata$factor.3.group,
                                  type1 = "mes_low",
                                  type2 = "mes_high",
                                  sig_scores = tcga_mofa_metadata[, comp_ids],
                                  path = file.path(out_path, "prrx1_and_factor_3_group/"),
                                  category_id = "factor_3_group",
                                  adjust_p_vals = FALSE,
                                  alpha = alpha,
                                  y_lab_vals = y_lab_vals,
                                  summary_list = summary_list)

# associations with immune signatures
summary_list <- list(tgfb = "TGFb_Family_Member")
y_lab_vals <- imm_n_map[colnames(imm_n_sig_scores_tcga_mofa)]
analyse_signature_by_focused_type(type = tcga_mofa_metadata$refined.type,
                                  type1 = "mes_low",
                                  type2 = "mes_high",
                                  sig_scores = imm_n_sig_scores_tcga_mofa,
                                  path = file.path(out_path, "immune_signatures_and_refined_types/"),
                                  category_id = "refined_type",
                                  alpha = alpha)

analyse_signature_by_focused_type(type = tcga_mofa_metadata$factor.3.group,
                                  type1 = "mes_low",
                                  type2 = "mes_high",
                                  sig_scores = imm_n_sig_scores_tcga_mofa,
                                  path = file.path(out_path, "immune_signatures_and_factor_3_group/"),
                                  category_id = "factor_3_group",
                                  alpha = alpha,
                                  y_lab_vals = y_lab_vals,
                                  summary_list = summary_list)

################################################################################
# analyse factor survival
################################################################################

surv_data <- cbind(tcga_mofa_factors,
                   tcga_mofa_factor_quantiles,
                   tcga_mofa_metadata,
                   tcga_mofa_patient_metadata)

# model survival associations with continuous factor scores
sp <- create_path(file.path(out_path, "mofa_factors_survival"))
var_ids <- colnames(tcga_mofa_factors)
analyse_survival("os_status",  sl["os"],  "os",  "os_months",  var_ids, surv_data, sp,
                 cov_ids1 = cov_ids1, samples_id = "mofa")
analyse_survival("dss_status", sl["dss"], "dss", "dss_months", var_ids, surv_data, sp,
                 cov_ids1 = cov_ids1, samples_id = "mofa")
analyse_survival("pfs_status", sl["pfs"], "pfs", "pfs_months", var_ids, surv_data, sp,
                 cov_ids1 = cov_ids1, samples_id = "mofa")

# model survival associations with factor 3 groups
sp <- create_path(file.path(out_path, "factor_3_group_survival"))
var_ids <- "factor.3.group"
analyse_survival("os_status",  sl["os"],  "os",  "os_months",  var_ids, surv_data, sp,
                 cov_ids1 = cov_ids1, samples_id = "mofa",
                 km_colors = refined_type_colors)
analyse_survival("dss_status", sl["dss"], "dss", "dss_months", var_ids, surv_data, sp,
                 cov_ids1 = cov_ids1, samples_id = "mofa",
                 km_colors = refined_type_colors)
analyse_survival("pfs_status", sl["pfs"], "pfs", "pfs_months", var_ids, surv_data, sp,
                 cov_ids1 = cov_ids1, samples_id = "mofa",
                 km_colors = refined_type_colors)

# These analyses do not converge
# model survival associations with refined types
#sp <- create_path(file.path(out_path, "refined_types_survival"))
#var_ids <- "refined.type"
#analyse_survival("os_status",  sl["os"],  "os",  "os_months",  var_ids, surv_data, sp,
#                 cov_ids1 = cov_ids1, cov_ids2 = cov_ids2, samples_id = "mofa")
#analyse_survival("dss_status", sl["dss"], "dss", "dss_months", var_ids, surv_data, sp,
#                 cov_ids1 = cov_ids1, cov_ids2 = cov_ids2, samples_id = "mofa")
#analyse_survival("pfs_status", sl["pfs"], "pfs", "pfs_months", var_ids, surv_data, sp,
#                 cov_ids1 = cov_ids1, cov_ids2 = cov_ids2, samples_id = "mofa")

# model survival associations with mofa types
#sp <- create_path(file.path(out_path, "mofa_types_survival"))
#var_ids <- "mofa.type"
#analyse_survival("os_status",  sl["os"],  "os",  "os_months",  var_ids, surv_data, sp,
#                 cov_ids1 = cov_ids1, cov_ids2 = cov_ids2, samples_id = "mofa")
#analyse_survival("dss_status", sl["dss"], "dss", "dss_months", var_ids, surv_data, sp,
#                 cov_ids1 = cov_ids1, cov_ids2 = cov_ids2, samples_id = "mofa")
#analyse_survival("pfs_status", sl["pfs"], "pfs", "pfs_months", var_ids, surv_data, sp,
#                 cov_ids1 = cov_ids1, cov_ids2 = cov_ids2, samples_id = "mofa")

################################################################################
# analyse factor-immune correlations
################################################################################

var_ids <- c(colnames(tcga_mofa_factors))

for(i in 1:length(var_ids)){
  trait_id <- var_ids[i]
  var_id <- tolower(trait_id)
  analyse_signature_by_trait(tcga_mofa_factors[,trait_id], 
                             imm_n_sig_scores_tcga_mofa, 
                             paste0(out_path, "/immune_signatures_and_mofa_factors/", var_id, "/"), 
                             trait_id, 
                             alpha=0.05)
}