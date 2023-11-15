
setwd("../")

################################################################################
# load library
################################################################################

library(ggplot2)
library(ggpubr)
source("../utilities/utils.R")
source("../utilities/signatures.R")
source("../utilities/mofa2_factor_calls.R")

################################################################################
# local functions
################################################################################

get_p_lab <- function(f.gg, y.col){
  stopifnot((!"y" %in% colnames(f.gg)) & is.data.frame(f.gg))
  f.gg$y <- f.gg[[y.col]]
  res <- lm(y ~ sh.and.dox, f.gg)
  p_lab <- paste0("P=", round(summary(res)$coefficients[2,4], 3))
  return(p_lab)
}

################################################################################
# load data
################################################################################

obj.expr <- load("R_Data/RNA_PRRX1_Data.RData")
obj.cell <- load("R_Data/RNA_Data.RData")
obj.sigs <- load("R_Data/Gene_Sigs.RData")
obj.mofa <- load(file = "R_Data/MOFA2_TopN_Heatmap_Feature_Filtering_Seed20200225_RNA_Weights.RData")
obj.meta <- load("R_Data/Heatmap_Metadata.RData")

rna_prrx1_metadata$sh <- sapply(rna_prrx1_metadata$shrna, function(x){
  if(x == "no_sh"){
    "no_sh"
  } else {
    "sh"
  }
})
rna_prrx1_metadata$sh.and.dox <- factor(rna_prrx1_metadata$sh == "sh" & rna_prrx1_metadata$condition == "DOX",
                                        levels = c("FALSE", "TRUE"))

prrx1_cl_rpkm_norm_uf <- cbind(rna_rpkm_norm_uf, 
                               rna_prrx1_rpkm_norm_uf[rownames(rna_rpkm_norm_uf),])

################################################################################
# plot signature changes
################################################################################

prrx1_sig_data <- sweep(rna_prrx1_rpkm_norm_uf, 1, rowMeans(rna_prrx1_rpkm_norm_uf))

tnbc_gsl_prrx1 <- uniform_subset_gsl(tnbc_gsl, rownames(prrx1_sig_data))

prrx1_sig_stats <- easy_signature_score(prrx1_sig_data, tnbc_gsl_prrx1)

# Hs578
gg <- cbind(prrx1_sig_stats, rna_prrx1_metadata)[1:8,]

pdf("analysis/rna_seq_prrx1_kd/signatures/hs578_prrx1_kd_bas_sig.pdf")
ggpaired(gg, x = "condition", y = "bas", id = "shrna", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "shrna", repel = TRUE,
         ylab = "Basal ignature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$bas), label = get_p_lab(gg, "bas"), hjust="left", size = 6) + 
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_kd/signatures/hs578_prrx1_kd_lum_sig.pdf")
ggpaired(gg, x = "condition", y = "lum", id = "shrna", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "shrna", repel = TRUE,
         ylab = "Luminal signature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$lum), label = get_p_lab(gg, "lum"), hjust="left", size = 6) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_kd/signatures/hs578_prrx1_kd_mes_sig.pdf")
ggpaired(gg, x = "condition", y = "mes", id = "shrna", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "shrna", repel = TRUE,
         ylab = "Mesenchymal signature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$mes), label = get_p_lab(gg, "mes"), hjust="left", size = 6) +
  my_theme()
dev.off()

# TT642
gg <- cbind(prrx1_sig_stats, rna_prrx1_metadata)[9:16,]

pdf("analysis/rna_seq_prrx1_kd/signatures/tt642_prrx1_kd_bas_sig.pdf")
ggpaired(gg, x = "condition", y = "bas", id = "shrna", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "shrna", repel = TRUE,
         ylab = "Basal ignature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$bas), label = get_p_lab(gg, "bas"), hjust="left", size = 6) + 
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_kd/signatures/tt642_prrx1_kd_lum_sig.pdf")
ggpaired(gg, x = "condition", y = "lum", id = "shrna", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "shrna", repel = TRUE,
         ylab = "Luminal signature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$lum), label = get_p_lab(gg, "lum"), hjust="left", size = 6) + 
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_kd/signatures/tt642_prrx1_kd_mes_sig.pdf")
ggpaired(gg, x = "condition", y = "mes", id = "shrna", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "shrna", repel = TRUE,
         ylab = "Mesenchymal signature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$mes), label = get_p_lab(gg, "mes"), hjust="left", size = 6) +
  my_theme()
dev.off()

################################################################################
# control - plot PRRX1 expression
################################################################################

gg <- cbind(PRRX1 = as.numeric(prrx1_sig_data["PRRX1",]), rna_prrx1_metadata)[1:8,]

pdf("analysis/rna_seq_prrx1_kd/signatures/hs578_prrx1_kd_prrx1_ctl.pdf")
ggpaired(gg, x = "condition", y = "PRRX1", id = "shrna", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "shrna", repel = TRUE,
         ylab = "PRRX1 Expression", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$PRRX1), label = get_p_lab(gg, "PRRX1"), hjust="left", size = 6) +
  my_theme()
dev.off()

gg <- cbind(PRRX1 = as.numeric(prrx1_sig_data["PRRX1",]), rna_prrx1_metadata)[9:16,]

pdf("analysis/rna_seq_prrx1_kd/signatures/tt642_prrx1_kd_prrx1_ctl.pdf")
ggpaired(gg, x = "condition", y = "PRRX1", id = "shrna", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "shrna", repel = TRUE,
         ylab = "PRRX1 Expression", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$PRRX1), label = get_p_lab(gg, "PRRX1"), hjust="left", size = 6) +
  my_theme()
dev.off()

################################################################################
# MOFA
################################################################################

stopifnot(all(sort(rownames(rna_rpkm_norm_uf)) == sort(rownames(rna_prrx1_rpkm_norm_uf))))
stopifnot(length(intersect(colnames(rna_rpkm_norm_uf), colnames(rna_prrx1_rpkm_norm_uf))) == 0)

keep <- rownames(rna_prrx1_rpkm_norm_uf) %in% rownames(rna_weights)
data_list <- list(rna = as.matrix(rna_prrx1_rpkm_norm_uf))
weight_list <- list(rna = rna_weights)
call_data <- prepare_mofa_calls(data_list, weight_list)
mofa_scores <- call_mofa_factor_scores(call_data)

gg <- cbind(mofa_scores, rna_prrx1_metadata)[1:8,]

pdf("analysis/rna_seq_prrx1_kd/signatures/hs578_prrx1_kd_f2_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor2", id = "shrna", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "shrna", repel = TRUE,
         ylab = "Factor 2", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor2), label = get_p_lab(gg, "Factor2"), hjust="left", size = 6) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_kd/signatures/hs578_prrx1_kd_f3_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor3", id = "shrna", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "shrna", repel = TRUE,
         ylab = "Factor 3", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor3), label = get_p_lab(gg, "Factor3"), hjust="left", size = 6) + 
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_kd/signatures/hs578_prrx1_kd_f6_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor6", id = "shrna", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "shrna", repel = TRUE,
         ylab = "Factor 6", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor6), label = get_p_lab(gg, "Factor6"), hjust="left", size = 6) +
  my_theme()
dev.off()

gg <- cbind(mofa_scores, rna_prrx1_metadata)[9:16,]

pdf("analysis/rna_seq_prrx1_kd/signatures/tt642_prrx1_kd_f2_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor2", id = "shrna", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "shrna", repel = TRUE,
         ylab = "Factor 2", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor2), label = get_p_lab(gg, "Factor2"), hjust="left", size = 6) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_kd/signatures/tt642_prrx1_kd_f3_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor3", id = "shrna", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "shrna", repel = TRUE,
         ylab = "Factor 3", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor3), label = get_p_lab(gg, "Factor3"), hjust="left", size = 6) + 
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_kd/signatures/tt642_prrx1_kd_f6_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor6", id = "shrna", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "shrna", repel = TRUE,
         ylab = "Factor 6", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor6), label = get_p_lab(gg, "Factor6"), hjust="left", size = 6) +
  my_theme()
dev.off()

keep <- rownames(prrx1_cl_rpkm_norm_uf) %in% rownames(rna_weights)
data_list <- list(rna = as.matrix(prrx1_cl_rpkm_norm_uf[keep, ]))
weight_list <- list(rna = rna_weights)
call_data <- prepare_mofa_calls(data_list, weight_list)
mofa_scores <- call_mofa_factor_scores(call_data)

gg <- data.frame(mofa_scores)
gg$exp <- c(rep("cl", nrow(rna_metadata)), rna_prrx1_metadata$sh)
gg$type <- c(rna_metadata$tnbc.subtype, rep(c("mesenchymal", "rhabdoid"), each = 8))
gg$condition <- c(rep("no_DOX", nrow(rna_metadata)), rna_prrx1_metadata$condition)

pdf("analysis/rna_seq_prrx1_kd/signatures/prrx1_kd_mofa_f2_f3_projection.pdf")
ggplot(gg, aes(x = Factor2, y = Factor3, color = type, shape = exp, fill = condition)) +
  geom_point() +
  scale_color_manual(values = c(line_colours$type, "rhabdoid" = "black")) +
  scale_shape_manual(values = c("cl" = 17, "no_sh" = 24, "sh" = 21)) +
  scale_fill_manual(values = c("no_DOX" = "white", "DOX" = "grey")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_kd/signatures/prrx1_kd_mofa_f2_f6_projection.pdf")
ggplot(gg, aes(x = Factor2, y = Factor6, color = type, shape = exp, fill = condition)) +
  geom_point() +
  scale_color_manual(values = c(line_colours$type, "rhabdoid" = "black")) +
  scale_shape_manual( values = c("cl" = 17, "no_sh" = 24, "sh" = 21)) +
  scale_fill_manual(values = c("no_DOX" = "white", "DOX" = "grey")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_kd/signatures/prrx1_kd_mofa_f3_f6_projection.pdf")
ggplot(gg, aes(x = Factor3, y = Factor6, color = type, shape = exp, fill = condition)) +
  geom_point() +
  scale_color_manual(values = c(line_colours$type, "rhabdoid" = "black")) +
  scale_shape_manual( values = c("cl" = 17, "no_sh" = 24, "sh" = 21)) +
  scale_fill_manual(values = c("no_DOX" = "white", "DOX" = "grey")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  my_theme()
dev.off()

################################################################################
# PRRX1 and FOXA1 control
################################################################################

gg <- data.frame(FOXA1 = as.numeric(prrx1_cl_rpkm_norm_uf["FOXA1",]),
                 PRRX1 = as.numeric(prrx1_cl_rpkm_norm_uf["PRRX1",]))

gg$exp <- c(rep("cl", nrow(rna_metadata)), rna_prrx1_metadata$sh)
gg$type <- c(rna_metadata$tnbc.subtype, rep(c("mesenchymal", "rhabdoid"), each = 8))
gg$condition <- c(rep("no_DOX", nrow(rna_metadata)), rna_prrx1_metadata$condition)

pdf("analysis/rna_seq_prrx1_kd/signatures/hs578_prrx1_kd_prrx1_foxa1_ctl.pdf")
ggplot(gg, aes(x = FOXA1, y = PRRX1, color = type, shape = exp, fill = condition)) +
  geom_point() +
  scale_color_manual(values = c(line_colours$type, "rhabdoid" = "black")) +
  scale_shape_manual(values = c("cl" = 17, "no_sh" = 24, "sh" = 21)) +
  scale_fill_manual(values = c("no_DOX" = "white", "DOX" = "grey")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  my_theme()
dev.off()
