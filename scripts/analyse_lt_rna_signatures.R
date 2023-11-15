
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

obj.expr <- load("R_Data/RNA_PRRX1_LT_Data.RData")
obj.cell <- load("R_Data/RNA_Data.RData")
obj.sigs <- load("R_Data/Gene_Sigs.RData")
obj.mofa <- load(file = "R_Data/MOFA2_TopN_Heatmap_Feature_Filtering_Seed20200225_RNA_Weights.RData")
obj.meta <- load("R_Data/Heatmap_Metadata.RData")

rna_lt_metadata$id <- paste0(rna_lt_metadata$shrna, "_", rna_lt_metadata$treatment.length)


rna_lt_metadata$sh <- sapply(rna_lt_metadata$shrna, function(x){
  if(x == "no_sh"){
    "no_sh"
  } else {
    "sh"
  }
})
rna_lt_metadata$sh.and.dox <- factor(rna_lt_metadata$sh == "sh" & rna_lt_metadata$condition == "DOX",
                                     levels = c("FALSE", "TRUE"))

lt_cl_rpkm_norm_uf <- cbind(rna_rpkm_norm_uf, 
                            rna_lt_rpkm_norm_uf[rownames(rna_rpkm_norm_uf),])

################################################################################
# plot signature changes
################################################################################

lt_sig_data <- sweep(rna_lt_rpkm_norm_uf, 1, rowMeans(rna_lt_rpkm_norm_uf))

tnbc_gsl_lt <- uniform_subset_gsl(tnbc_gsl, rownames(lt_sig_data))

lt_sig_stats <- easy_signature_score(lt_sig_data, tnbc_gsl_lt)

# Hs578 4 weeks
gg <- subset(cbind(lt_sig_stats, rna_lt_metadata), parental.cell.line == "Hs578" & treatment.length == "4 weeks")

pdf("analysis/rna_seq_prrx1_lt/signatures/hs578_prrx1_4w_bas_sig.pdf")
ggpaired(gg, x = "condition", y = "bas", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Basal signature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$bas), label = get_p_lab(gg, "bas"), hjust="left", size = 6) + 
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/hs578_prrx1_4w_lum_sig.pdf")
ggpaired(gg, x = "condition", y = "lum", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Luminal signature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$lum), label = get_p_lab(gg, "lum"), hjust="left", size = 6) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/hs578_prrx1_4w_mes_sig.pdf")
ggpaired(gg, x = "condition", y = "mes", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Mesenchymal signature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$mes), label = get_p_lab(gg, "mes"), hjust="left", size = 6) +
  my_theme()
dev.off()

# Hs578 8 weeks
gg <- subset(cbind(lt_sig_stats, rna_lt_metadata), parental.cell.line == "Hs578" & treatment.length == "8 weeks")

pdf("analysis/rna_seq_prrx1_lt/signatures/hs578_prrx1_8w_bas_sig.pdf")
ggpaired(gg, x = "condition", y = "bas", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Basal signature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$bas), label = get_p_lab(gg, "bas"), hjust="left", size = 6) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/hs578_prrx1_8w_lum_sig.pdf")
ggpaired(gg, x = "condition", y = "lum", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Luminal signature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$lum), label = get_p_lab(gg, "lum"), hjust="left", size = 6) + 
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/hs578_prrx1_8w_mes_sig.pdf")
ggpaired(gg, x = "condition", y = "mes", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Mesenchymal signature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$mes), label = get_p_lab(gg, "mes"), hjust="left", size = 6) +
  my_theme()
dev.off()

# TT642 4 weeks
gg <- subset(cbind(lt_sig_stats, rna_lt_metadata), parental.cell.line == "TT642" & treatment.length == "4 weeks")

pdf("analysis/rna_seq_prrx1_lt/signatures/tt642_prrx1_4w_bas_sig.pdf")
ggpaired(gg, x = "condition", y = "bas", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Basal signature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$bas), label = get_p_lab(gg, "bas"), hjust="left", size = 6) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/tt642_prrx1_4w_lum_sig.pdf")
ggpaired(gg, x = "condition", y = "lum", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Luminal signature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$lum), label = get_p_lab(gg, "lum"), hjust="left", size = 6) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/tt642_prrx1_4w_mes_sig.pdf")
ggpaired(gg, x = "condition", y = "mes", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Mesenchymal ignature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$mes), label = get_p_lab(gg, "mes"), hjust="left", size = 6) +
  my_theme()
dev.off()

# TT642 8 weeks
gg <- subset(cbind(lt_sig_stats, rna_lt_metadata), parental.cell.line == "TT642" & treatment.length == "8 weeks")

pdf("analysis/rna_seq_prrx1_lt/signatures/tt642_prrx1_8w_bas_sig.pdf")
ggpaired(gg, x = "condition", y = "bas", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Basal signature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$bas), label = get_p_lab(gg, "bas"), hjust="left", size = 6) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/tt642_prrx1_8w_lum_sig.pdf")
ggpaired(gg, x = "condition", y = "lum", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Luminal signature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$lum), label = get_p_lab(gg, "lum"), hjust="left", size = 6) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/tt642_prrx1_8w_mes_sig.pdf")
ggpaired(gg, x = "condition", y = "mes", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Mesenchymal ignature", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$mes), label = get_p_lab(gg, "mes"), hjust="left", size = 6) +
  my_theme()
dev.off()

################################################################################
# control - plot PRRX1 expression
################################################################################

gg <- subset(cbind(PRRX1 = as.numeric(lt_sig_data["PRRX1",]), rna_lt_metadata), parental.cell.line == "Hs578" & treatment.length == "4 weeks")

pdf("analysis/rna_seq_prrx1_lt/signatures/hs578_prrx1_4w_prrx1_ctl.pdf")
ggpaired(gg, x = "condition", y = "PRRX1", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "PRRX1 Expression", font.label = list(size = 20)) + 
  annotate("text", 0.5, max(gg$PRRX1), label = get_p_lab(gg, "PRRX1"), hjust="left", size = 6) +
  my_theme()
dev.off()

gg <- subset(cbind(PRRX1 = as.numeric(lt_sig_data["PRRX1",]), rna_lt_metadata), parental.cell.line == "Hs578" & treatment.length == "8 weeks")

pdf("analysis/rna_seq_prrx1_lt/signatures/hs578_prrx1_8w_prrx1_ctl.pdf")
ggpaired(gg, x = "condition", y = "PRRX1", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "PRRX1 Expression", font.label = list(size = 20)) + 
  annotate("text", 0.5, max(gg$PRRX1), label = get_p_lab(gg, "PRRX1"), hjust="left", size = 6) +
  my_theme()
dev.off()

gg <- subset(cbind(PRRX1 = as.numeric(lt_sig_data["PRRX1",]), rna_lt_metadata), parental.cell.line == "TT642" & treatment.length == "4 weeks")

pdf("analysis/rna_seq_prrx1_lt/signatures/tt642_prrx1_4w_prrx1_ctl.pdf")
ggpaired(gg, x = "condition", y = "PRRX1", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "PRRX1 Expression", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$PRRX1), label = get_p_lab(gg, "PRRX1"), hjust="left", size = 6) +
  my_theme()
dev.off()

gg <- subset(cbind(PRRX1 = as.numeric(lt_sig_data["PRRX1",]), rna_lt_metadata), parental.cell.line == "TT642" & treatment.length == "8 weeks")

pdf("analysis/rna_seq_prrx1_lt/signatures/tt642_prrx1_8w_prrx1_ctl.pdf")
ggpaired(gg, x = "condition", y = "PRRX1", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "PRRX1 Expression", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$PRRX1), label = get_p_lab(gg, "PRRX1"), hjust="left", size = 6) +
  my_theme()
dev.off()

################################################################################
# MOFA projection
################################################################################

stopifnot(all(sort(rownames(rna_rpkm_norm_uf)) == sort(rownames(rna_lt_rpkm_norm_uf))))
stopifnot(length(intersect(colnames(rna_rpkm_norm_uf), colnames(rna_lt_rpkm_norm_uf))) == 0)

keep <- rownames(rna_lt_rpkm_norm_uf) %in% rownames(rna_weights)
data_list <- list(rna = as.matrix(rna_lt_rpkm_norm_uf))
weight_list <- list(rna = rna_weights)
call_data <- prepare_mofa_calls(data_list, weight_list)
mofa_scores <- call_mofa_factor_scores(call_data)

# Hs578 4w
gg <- subset(cbind(mofa_scores, rna_lt_metadata), parental.cell.line == "Hs578" & treatment.length == "4 weeks")

pdf("analysis/rna_seq_prrx1_lt/signatures/hs578_prrx1_4w_f2_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor2", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Factor 2", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor2), label = get_p_lab(gg, "Factor2"), hjust="left", size = 6) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/hs578_prrx1_4w_f3_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor3", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Factor 3", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor3), label = get_p_lab(gg, "Factor3"), hjust="left", size = 6) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/hs578_prrx1_4w_f6_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor6", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Factor 6", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor6), label = get_p_lab(gg, "Factor6"), hjust="left", size = 6) +
  my_theme()
dev.off()

# Hs578 8w
gg <- subset(cbind(mofa_scores, rna_lt_metadata), parental.cell.line == "Hs578" & treatment.length == "8 weeks")

pdf("analysis/rna_seq_prrx1_lt/signatures/hs578_prrx1_8w_f2_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor2", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Factor 2", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor2), label = get_p_lab(gg, "Factor2"), hjust="left", size = 6) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/hs578_prrx1_8w_f3_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor3", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Factor 3", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor3), label = get_p_lab(gg, "Factor3"), hjust="left", size = 6) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/hs578_prrx1_8w_f6_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor6", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Factor 6", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor6), label = get_p_lab(gg, "Factor6"), hjust="left", size = 6) +
  my_theme()
dev.off()

# TT642 4w
gg <- subset(cbind(mofa_scores, rna_lt_metadata), parental.cell.line == "TT642" & treatment.length == "4 weeks")

pdf("analysis/rna_seq_prrx1_lt/signatures/tt642_prrx1_4w_f2_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor2", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Factor 2", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor2), label = get_p_lab(gg, "Factor2"), hjust="left", size = 6) + 
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/tt642_prrx1_4w_f3_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor3", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Factor 3", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor3), label = get_p_lab(gg, "Factor3"), hjust="left", size = 6) + 
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/tt642_prrx1_4w_f6_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor6", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Factor 6", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor6), label = get_p_lab(gg, "Factor6"), hjust="left", size = 6) + 
  my_theme()
dev.off()

# TT642 8w
gg <- subset(cbind(mofa_scores, rna_lt_metadata), parental.cell.line == "TT642" & treatment.length == "8 weeks")

pdf("analysis/rna_seq_prrx1_lt/signatures/tt642_prrx1_8w_f2_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor2", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Factor 2", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor2), label = get_p_lab(gg, "Factor2"), hjust="left", size = 6) + 
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/tt642_prrx1_8w_f3_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor3", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Factor 3", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor3), label = get_p_lab(gg, "Factor3"), hjust="left", size = 6) + 
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/tt642_prrx1_8w_f6_sig.pdf")
ggpaired(gg, x = "condition", y = "Factor6", id = "id", 
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "jco", label = "id", repel = TRUE,
         ylab = "Factor 6", font.label = list(size = 20)) +
  annotate("text", 0.5, max(gg$Factor6), label = get_p_lab(gg, "Factor6"), hjust="left", size = 6) + 
  my_theme()
dev.off()

keep <- rownames(lt_cl_rpkm_norm_uf) %in% rownames(rna_weights)
data_list <- list(rna = as.matrix(lt_cl_rpkm_norm_uf[keep, ]))
weight_list <- list(rna = rna_weights)
call_data <- prepare_mofa_calls(data_list, weight_list)
mofa_scores <- call_mofa_factor_scores(call_data)

gg <- data.frame(mofa_scores)
gg$exp <- c(rep("cl", nrow(rna_metadata)), rna_lt_metadata$sh)
gg$type <- c(rna_metadata$tnbc.subtype, rep("mesenchymal", 10), rep("rhabdoid", 12))
gg$condition <- c(rep("no_DOX", nrow(rna_metadata)), rna_lt_metadata$condition)
ii_1 <- c(1:nrow(rna_metadata), nrow(rna_metadata) + which(rna_lt_metadata$treatment.length == "4 weeks"))
ii_2 <- c(1:nrow(rna_metadata), nrow(rna_metadata) + which(rna_lt_metadata$treatment.length == "8 weeks"))

# 4 weeks
pdf("analysis/rna_seq_prrx1_lt/signatures/prrx1_4w_mofa_f2_f3_projection.pdf")
ggplot(gg[ii_1,], aes(x = Factor2, y = Factor3, color = type, shape = exp, fill = condition)) +
  geom_point() +
  scale_color_manual(values = c(line_colours$type, "rhabdoid" = "black")) +
  scale_shape_manual( values = c("cl" = 17, "no_sh" = 24, "sh" = 21)) +
  scale_fill_manual(values = c("no_DOX" = "white", "DOX" = "grey")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/prrx1_4w_mofa_f2_f6_projection.pdf")
ggplot(gg[ii_1,], aes(x = Factor2, y = Factor6, color = type, shape = exp, fill = condition)) +
  geom_point() +
  scale_color_manual(values = c(line_colours$type, "rhabdoid" = "black")) +
  scale_shape_manual( values = c("cl" = 17, "no_sh" = 24, "sh" = 21)) +
  scale_fill_manual(values = c("no_DOX" = "white", "DOX" = "grey")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/prrx1_4w_mofa_f3_f6_projection.pdf")
ggplot(gg[ii_1,], aes(x = Factor3, y = Factor6, color = type, shape = exp, fill = condition)) +
  geom_point() +
  scale_color_manual(values = c(line_colours$type, "rhabdoid" = "black")) +
  scale_shape_manual( values = c("cl" = 17, "no_sh" = 24, "sh" = 21)) +
  scale_fill_manual(values = c("no_DOX" = "white", "DOX" = "grey")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  my_theme()
dev.off()

# 8 weeks
pdf("analysis/rna_seq_prrx1_lt/signatures/prrx1_8w_mofa_f2_f3_projection.pdf")
ggplot(gg[ii_2,], aes(x = Factor2, y = Factor3, color = type, shape = exp, fill = condition)) +
  geom_point() +
  scale_color_manual(values = c(line_colours$type, "rhabdoid" = "black")) +
  scale_shape_manual( values = c("cl" = 17, "no_sh" = 24, "sh" = 21)) +
  scale_fill_manual(values = c("no_DOX" = "white", "DOX" = "grey")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/prrx1_8w_mofa_f2_f6_projection.pdf")
ggplot(gg[ii_2,], aes(x = Factor2, y = Factor6, color = type, shape = exp, fill = condition)) +
  geom_point() +
  scale_color_manual(values = c(line_colours$type, "rhabdoid" = "black")) +
  scale_shape_manual( values = c("cl" = 17, "no_sh" = 24, "sh" = 21)) +
  scale_fill_manual(values = c("no_DOX" = "white", "DOX" = "grey")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/prrx1_8w_mofa_f3_f6_projection.pdf")
ggplot(gg[ii_2,], aes(x = Factor3, y = Factor6, color = type, shape = exp, fill = condition)) +
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

gg <- data.frame(FOXA1 = as.numeric(lt_cl_rpkm_norm_uf["FOXA1",]),
                 PRRX1 = as.numeric(lt_cl_rpkm_norm_uf["PRRX1",]))

gg$exp <- c(rep("cl", nrow(rna_metadata)), rna_lt_metadata$sh)
gg$type <- c(rna_metadata$tnbc.subtype, rep("mesenchymal", 10), rep("rhabdoid", 12))
gg$condition <- c(rep("no_DOX", nrow(rna_metadata)), rna_lt_metadata$condition)
ii_1 <- c(1:nrow(rna_metadata), nrow(rna_metadata) + which(rna_lt_metadata$treatment.length == "4 weeks"))
ii_2 <- c(1:nrow(rna_metadata), nrow(rna_metadata) + which(rna_lt_metadata$treatment.length == "8 weeks"))

pdf("analysis/rna_seq_prrx1_lt/signatures/hs578_prrx1_4w_prrx1_foxa1_ctl.pdf")
ggplot(gg[ii_1,], aes(x = FOXA1, y = PRRX1, color = type, shape = exp, fill = condition)) +
  geom_point() +
  scale_color_manual(values = c(line_colours$type, "rhabdoid" = "black")) +
  scale_shape_manual( values = c("cl" = 17, "no_sh" = 24, "sh" = 21)) +
  scale_fill_manual(values = c("no_DOX" = "white", "DOX" = "grey")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  my_theme()
dev.off()

pdf("analysis/rna_seq_prrx1_lt/signatures/hs578_prrx1_8w_prrx1_foxa1_ctl.pdf")
ggplot(gg[ii_2,], aes(x = FOXA1, y = PRRX1, color = type, shape = exp, fill = condition)) +
  geom_point() +
  scale_color_manual(values = c(line_colours$type, "rhabdoid" = "black")) +
  scale_shape_manual( values = c("cl" = 17, "no_sh" = 24, "sh" = 21)) +
  scale_fill_manual(values = c("no_DOX" = "white", "DOX" = "grey")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  my_theme()
dev.off()
