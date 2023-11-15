source("../utilities/utils.R")
source("../utilities/sampling.R")
source("../utilities/signatures.R")
source("../utilities/mofa2_factor_calls.R")
source("../src/hexplot/hexplot.R")
library(Seurat)
library(reshape2)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(SingleCellExperiment)
library(foreach)
library(doMC)
library(mclust)
registerDoMC(22)

setwd("../")

# Load data ---------------------------------------------------------------------------------------------------------------------------------------

load(file = "R_Data/SC_BC_DR.RData")
load(file = "R_Data/SC_Metadata.RData")
load(file = "R_Data/Heatmap_Metadata.RData")
load(file = "R_Data/Gene_Sigs.RData")
load(file = "R_Data/RNA_Data.RData")
load(file = "R_Data/RNA_Metadata.RData")
load(file = "R_Data/MOFA2_TopN_Heatmap_Feature_Filtering_Seed20200225_RNA_Weights.RData")

line_sn_colours <- line_colours
names(line_sn_colours$type) <- c("bas", "mes", "lum")
hex_colours <- c("none" = "#808080", line_sn_colours$type[c("bas", "lum", "mes")])
hex_colours[5] <- easy.midColor(hex_colours[c(2,3)])
hex_colours[6] <- easy.midColor(hex_colours[c(2,4)])
hex_colours[7] <- easy.midColor(hex_colours[c(3,4)])
hex_colours[8] <- easy.midColor(c(hex_colours[2], easy.midColor(hex_colours[c(3,4)])))
names(hex_colours)[5:8] <- c("bas_and_lum", "bas_and_mes", "lum_and_mes", "all")

# TNBC Signatures ---------------------------------------------------------------------------------------------------------------------------------------

sc$parental.tnbc.type <- line_annotation[sc$orig.ident, "type"]

pdf("analysis/scrna_seq/signatures/parental_tnbc_type.pdf")
DimPlot(sc, group.by = "parental.tnbc.type", cols = line_colours$type)
dev.off()

pdf("analysis/scrna_seq/signatures/vimentin.pdf")
FeaturePlot(sc, "VIM", label = TRUE)
dev.off()

pdf("analysis/scrna_seq/other_analyses/igfbp7.pdf")
FeaturePlot(sc, "IGFBP7")
dev.off()

# Project onto MOFA tnbc type-linked factors
if(file.exists("R_Data/SC_MOFA_Scores.RData")){
  load("R_Data/SC_MOFA_Scores.RData")
} else{
  keep <- rownames(sc) %in% rownames(rna_weights)
  data_list <- list(rna = sc[["lnbc"]]@data[keep, ])
  weight_list <- list(rna = rna_weights)
  call_data <- prepare_mofa_calls(data_list, weight_list)
  mofa_scores <- call_mofa_factor_scores(call_data)
  orig_idents <- sc$orig.ident
  save(orig_idents, mofa_scores, file = "R_Data/SC_MOFA_Scores.RData")
}

# MOFA ITH analysis
cl <- factor(sc@meta.data$orig.ident)
wss <- apply(mofa_scores, 2, function(x){
  sum(tapply(x, cl, var) * (table(cl) - 1))
})
tss <- apply(mofa_scores, 2, var) * (nrow(mofa_scores) - 1)
mofa_het <- data.frame(wss = wss, tss = tss)
mofa_het$score <- mofa_het$wss / mofa_het$tss
mofa_het$factor <- colnames(mofa_scores)

pp <- mofa_het
pp$factor <- factor(pp$factor, levels = pp$factor[order(pp$score)])
pdf("analysis/scrna_seq/mofa_projection/within_ss.pdf")
ggplot(pp, aes(x = factor, y = score, fill = factor)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  ylab("Intra-cell line variance proportion") + 
  my_theme()
dev.off()

# Assess evidence for MOFA score clusters
if(file.exists("R_Data/SC_MOFA_LL.RData")){
  load("R_Data/SC_MOFA_LL.RData")
} else {
  ll <- list()
  for(i in 1:nrow(scrna_sample_table)){
    sample <- as.character(scrna_sample_table$sample3[i])
    print(sample)
    mm <- mofa_scores[sc$sample3 == sample, c(2, 3, 6)]
    ll[[sample]] <- sapply(1:5, function(x) Mclust(mm, x)$loglik)
  }
  save(ll, file = "R_Data/SC_MOFA_LL.RData")  
}

ll2 <- sapply(ll, function(x) x[2:5]-x[1:4])
p <- melt(cbind(as.data.frame(ll2), k = 2:5), id.var = "k")
colnames(p) <- c("Clusters", "Sample", "DeltaLogLik") 
p$Sample <- factor(p$Sample, levels = sort(scrna_sample_table$sample3))
p$Label <- as.character(p$Sample)
p$Label[which(p$Clusters != 2 | !p$Sample %in% c("HDQP1", "SUM185", "EMG3"))] = ""
pdf("analysis/scrna_seq/mofa_projection/clustering/delta_loglik_vs_last.pdf")
ggplot(p, aes(x = Clusters, y = DeltaLogLik, color = Sample)) + 
  geom_point() +
  geom_text_repel(aes(label = Label), max.overlaps = Inf, show.legend = FALSE) +
  ylab("Change in log likelihood") +
  my_theme() +
  theme(legend.position="top", legend.title = element_blank())
dev.off()

x <- as.matrix(mofa_scores)
colnames(x) <- paste0("MOFA_", 1:ncol(x))
sc[["mofa"]] <- CreateDimReducObject(embeddings = x, key = "MOFA_", assay = "lnbc")

for(i in 1:nrow(scrna_sample_table)){
  sample <- as.character(scrna_sample_table$sample3[i])
  print(sample)
  sc@meta.data[[sample]] <- sc$sample3 == sample
  
  pdf_prefix <- paste0("analysis/scrna_seq/mofa_projection/by_cell_line/sample_", sample)
  pdf(paste0(pdf_prefix, "_mofa_factor_2_and_3_projection.pdf"))
  p <- DimPlot(sc, reduction = "mofa", group.by = sample, dims = c(2, 3))
  print(p)
  dev.off()
  pdf(paste0(pdf_prefix, "_mofa_factor_2_and_6_projection.pdf"))
  p <- DimPlot(sc, reduction = "mofa", group.by = sample, dims = c(2, 6))
  print(p)
  dev.off()
  pdf(paste0(pdf_prefix, "_mofa_factor_3_and_6_projection.pdf"))
  p <- DimPlot(sc, reduction = "mofa", group.by = sample, dims = c(3, 6))
  print(p)
  dev.off()
  
  x <- sc[,sc@meta.data[[sample]]]
  pdf_prefix <- paste0("analysis/scrna_seq/mofa_projection/by_cell_line/focused_sample_", sample)
  pdf(paste0(pdf_prefix, "_mofa_factor_2_and_3_projection.pdf"))
  p <- DimPlot(x, reduction = "mofa", group.by = sample, dims = c(2, 3))
  print(p)
  dev.off()
  pdf(paste0(pdf_prefix, "_mofa_factor_2_and_6_projection.pdf"))
  p <- DimPlot(x, reduction = "mofa", group.by = sample, dims = c(2, 6))
  print(p)
  dev.off()
  pdf(paste0(pdf_prefix, "_mofa_factor_3_and_6_projection.pdf"))
  p <- DimPlot(x, reduction = "mofa", group.by = sample, dims = c(3, 6))
  print(p)
  dev.off()
}

mofa_scores_df <- as.data.frame(mofa_scores)
mofa_scores_df$sample <- sc$sample3
mofa_scores_agg <- aggregate(. ~ sample, mofa_scores_df, mean)
mofa_scores_agg$tnbc.type <- scrna_sample_table$tnbc.subtype[match(mofa_scores_agg$sample, scrna_sample_table$sample3)]

pdf("analysis/scrna_seq/mofa_projection/mofa_factor_2_and_3_mean_pos.pdf")
ggplot(data = mofa_scores_agg, aes(x = Factor2, y = Factor3)) +
  geom_point(aes(color = tnbc.type)) +
  scale_color_manual(values = line_colours$type) +
  geom_text_repel(aes(label = sample)) +
  xlab("Factor 2 Score") +
  ylab("Factor 3 Score") +
  NoLegend() 
dev.off()

pdf("analysis/scrna_seq/mofa_projection/mofa_factor_2_and_6_mean_pos.pdf")
ggplot(data = mofa_scores_agg, aes(x = Factor2, y = Factor6)) +
  geom_point(aes(color = tnbc.type)) +
  scale_color_manual(values = line_colours$type) +
  geom_text_repel(aes(label = sample)) +
  xlab("Factor 2 Score") +
  ylab("Factor 6 Score") +
  NoLegend() 
dev.off()

pdf("analysis/scrna_seq/mofa_projection/mofa_factor_3_and_6_mean_pos.pdf")
ggplot(data = mofa_scores_agg, aes(x = Factor3, y = Factor6)) +
  geom_point(aes(color = tnbc.type)) +
  scale_color_manual(values = line_colours$type) +
  geom_text_repel(aes(label = sample)) +
  xlab("Factor 3 Score") +
  ylab("Factor 6 Score") +
  NoLegend() 
dev.off()

# Cell cycle ----------------------------------------------------------------------------------------------------

sc <- CellCycleScoring(sc, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

pdf("analysis/scrna_seq/signatures/cell_cycle_phase.pdf")
DimPlot(sc, group.by = "Phase")
dev.off()

phase_by_type <- table(sc$Phase, sc$parental.tnbc.type)
df <- melt(sweep(phase_by_type, 2, colSums(phase_by_type), "/"))
colnames(df) <- c("phase", "tnbc.type", "proportion")
pdf("analysis/scrna_seq/signatures/cell_cycle_by_type.pdf")
ggplot(df, aes(x = tnbc.type, y = proportion, fill = phase)) +
  geom_bar(stat = "identity", position = position_dodge())
dev.off()

phase_by_line <- table(sc$Phase, sc$sample3)
cc_prop <- apply(phase_by_line, 2, function(x) x / sum(x))
df <- melt(cc_prop)
colnames(df) <- c("phase", "sample", "proportion")
pdf("analysis/scrna_seq/signatures/cell_cycle_by_line.pdf")
ggplot(df, aes(x = sample, y = proportion, fill = phase)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


# Hexplots ----------------------------------------------------------------------------------------------------------------

# down-sample
n_samp <- min(table(sc$sample3))
ids <- BalancedSample(sc$sample3, n_samp, 20201029)

sc.ds <- sc[,ids]
tnbc_gsl_sub <- uniform_subset_gsl(tnbc_gsl, rownames(sc.ds))

samples <- unique(sc.ds$sample3)

exp_data <- as.matrix(sc.ds[["lnbc"]]@scale.data)

props <- list()
for(i in 1:length(samples)){
  sample <- samples[i]
  print(sample)
  
  ed <- exp_data[, which(sc.ds$sample3 == sample)]
  
  sig_score <- local.gsl.permutation(ed, tnbc_gsl_sub)
  
  pdf(paste0("analysis/scrna_seq/hexplot/scaled_", sample, "_hexplot.pdf"), 2.5, 2.5)
  easy.hexplot(sig_score$perm.pval, 
               no.text = FALSE, 
               cex = 0.5, 
               lwd = 0.2, 
               xcol = c("#808080", line_colours$type[c("basal", "luminal", "mesenchymal")]))
  dev.off()
  
  res <- prob2type(sig_score$perm.pval)
  res <- table(factor(res, levels = names(hex_colours)))
  props[[i]] <- res/sum(res)
}
names(props) <- samples
drop_types <- names(which(colMeans(do.call(rbind, props)) == 0))

props_long <- melt(props)
colnames(props_long) <- c("signature", "proportion", "sample")
props_long <- subset(props_long, !(signature %in% drop_types))
props_long$sample <- factor(props_long$sample, levels = levels(scrna_sample_table$sample3))
props_long$signature <- factor(props_long$signature, levels = names(hex_colours))

pdf("analysis/scrna_seq/hexplot/hexplot_summary_barplot.pdf")
ggplot(props_long, aes(x = sample, y = proportion, fill = signature)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = hex_colours) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Weighted version
exp_data <- as.matrix(sc.ds[["lnbc"]]@data)
exp_data <- sweep(exp_data, 1, rowMeans(exp_data))

weighted_props <- list()
for(i in 1:length(samples)){
  sample <- samples[i]
  print(sample)
  
  ed <- exp_data[, which(sc.ds$sample3 == sample)]
  
  sig_score <- local.gsl.permutation(ed, tnbc_gsl_sub)
  pdf(paste0("analysis/scrna_seq/hexplot/weighted_", sample, "_hexplot.pdf"), 2.5, 2.5)
  easy.hexplot(sig_score$perm.pval,
               no.text = FALSE,
               cex = 0.5,
               lwd = 0.2,
               xcol = c("#808080", line_colours$type[c("basal", "luminal", "mesenchymal")]))
  dev.off()
  
  res <- prob2type(sig_score$perm.pval)
  res <- table(factor(res, levels = names(hex_colours)))
  weighted_props[[i]] <- res/sum(res)
}
names(weighted_props) <- samples
wt_drop_types <- names(which(colMeans(do.call(rbind, weighted_props)) == 0))

wt_props_long <- melt(weighted_props)
colnames(wt_props_long) <- c("signature", "proportion", "sample")
wt_props_long <- subset(wt_props_long, !(signature %in% wt_drop_types))
wt_props_long$sample <- factor(wt_props_long$sample, levels = levels(scrna_sample_table$sample3))
wt_props_long$signature <- factor(wt_props_long$signature, levels = names(hex_colours))

pdf("analysis/scrna_seq/hexplot/weighted_hexplot_summary_barplot.pdf")
ggplot(wt_props_long, aes(x = sample, y = proportion, fill = signature)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = hex_colours) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

save(props, weighted_props, file = "R_Data/SC_Hex_Types.RData")
