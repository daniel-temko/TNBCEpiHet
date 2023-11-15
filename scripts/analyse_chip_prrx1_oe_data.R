# Code to analyse influence of PRRX1 over-expression on super-enhancer landscape in 
# basal and luminal cell lines
setwd("../")

##########################################################################################
## Load library
##########################################################################################

source("scripts/utilities/preprocessing.R")
source("../utilities/mofa2_factor_calls.R")
source("../utilities/utils.R")
library(pheatmap)
library(ggplot2)
library(ggrepel)

################################################################################
# local functions
################################################################################

up_tfs <- function(df){
  df[which((df$log2FoldChange > 0) & (df$dbd.tf)),]
}

up_ut_tfs <- function(df, tf_col){
  df[which((df$log2FoldChange > 0) & (df$ut.tf)),]
}

##########################################################################################
## Load data
##########################################################################################

obj.meta <- load("R_Data/Heatmap_Metadata.RData")
obj.mofa <- load(file = "R_Data/MOFA2_TopN_Heatmap_Feature_Filtering_Seed20200225_ChIP_Weights.RData")
obj.chip <- load("R_Data/ChIP_Data.RData")
obj.chip.oe <- load("R_Data/ChIP_PRRX1_OE_Data.RData")
cap_val <- 3
col_map <- list(type = line_colours$type2, 
                source = c("untreated" = "grey", "sh" = "orange"),
                mutant = c("wt" = "lightgrey", "dH3" = "violet"),
                time.point = c("short" = "purple", "long" = "black"),
                condition = c("no_dox" = "pink", "dox" = "yellow"))

bas_de <- read_csv("analysis/rna_seq/differential_genes/basal_de_lfc0_genes.csv")
mes_de <- read_csv("analysis/rna_seq/differential_genes/mesenchymal_de_lfc0_genes.csv")
lum_de <- read_csv("analysis/rna_seq/differential_genes/luminal_de_lfc0_genes.csv")

bas_dse <- read_csv("analysis/chip_seq/differential_regions/superenhancers_basal_diff_lfc0.csv")
mes_dse <- read_csv("analysis/chip_seq/differential_regions/superenhancers_mesenchymal_diff_lfc0.csv")
lum_dse <- read_csv("analysis/chip_seq/differential_regions/superenhancers_luminal_diff_lfc0.csv")

# Note: The logic here could be simplified, defining type tfs is not necessary
bas_tfs <- intersect(up_tfs(bas_de)$gene.symbol, up_tfs(bas_dse)$closest.gene)
mes_tfs <- intersect(up_tfs(mes_de)$gene.symbol, up_tfs(mes_dse)$closest.gene)
lum_tfs <- intersect(up_tfs(lum_de)$gene.symbol, up_tfs(lum_dse)$closest.gene)

bas_ses <- subset(up_tfs(bas_dse), closest.gene %in% bas_tfs)$region
mes_ses <- subset(up_tfs(mes_dse), closest.gene %in% mes_tfs)$region
lum_ses <- subset(up_tfs(lum_dse), closest.gene %in% lum_tfs)$region

bas_ut_tfs <- intersect(up_ut_tfs(bas_de)$gene.symbol, up_ut_tfs(bas_dse)$closest.gene)
mes_ut_tfs <- intersect(up_ut_tfs(mes_de)$gene.symbol, up_ut_tfs(mes_dse)$closest.gene)
lum_ut_tfs <- intersect(up_ut_tfs(lum_de)$gene.symbol, up_ut_tfs(lum_dse)$closest.gene)

bas_ut_ses <- subset(up_ut_tfs(bas_dse), closest.gene %in% bas_ut_tfs)$region
mes_ut_ses <- subset(up_ut_tfs(mes_dse), closest.gene %in% mes_ut_tfs)$region
lum_ut_ses <- subset(up_ut_tfs(lum_dse), closest.gene %in% lum_ut_tfs)$region

##########################################################################################
## Transform data
##########################################################################################

chip_oe_metadata$parental.type <- line_annotation[chip_oe_metadata$parental.cell.line, "type"]

# use informative sample names for OE data
chip_oe_metadata$sample.name <- paste0(chip_oe_metadata$parental.cell.line, "_", rep(1:16, 2))
rownames(chip_oe_metadata) <- chip_oe_metadata$sample.name
colnames(chip_oe_se_all_rpkm_norm) <- rownames(chip_oe_metadata)
colnames(chip_oe_se_all_rpkm_norm_uf) <- rownames(chip_oe_metadata)
colnames(chip_oe_se_all_rpkm) <- rownames(chip_oe_metadata)

stopifnot(all(rownames(chip_se_all_rpkm_norm_uf) == rownames(chip_oe_se_all_rpkm_norm_uf)))
comb_se_data <- cbind(chip_se_all_rpkm_norm_uf, chip_oe_se_all_rpkm_norm_uf)
comb_se_data_centered <- sweep(comb_se_data, 1, rowMeans(comb_se_data))

comb_metadata <- data.frame(row.names = c(rownames(chip_metadata), rownames(chip_oe_metadata)),
                            type = c(chip_metadata$tnbc.subtype, chip_oe_metadata$parental.type),
                            source = c(rep("untreated", nrow(chip_metadata)), rep("sh", nrow(chip_oe_metadata))),
                            mutant = c(rep(NA, nrow(chip_metadata)), chip_oe_metadata$mutant),
                            time.point = c(rep(NA, nrow(chip_metadata)), chip_oe_metadata$time.point),
                            condition = c(rep(NA, nrow(chip_metadata)), chip_oe_metadata$condition))
comb_metadata$sample <- rownames(comb_metadata)
comb_metadata$label <- sapply(1:nrow(comb_metadata), function(x){
  if(comb_metadata$source[x] == "untreated"){
    comb_metadata$sample[x]
  } else {
    ""
  }
})
comb_metadata$group <- sapply(1:nrow(comb_metadata), function(x){
  if(comb_metadata$source[x] == "untreated"){
    "untreated"
  } else{
    paste(comb_metadata$source[x], comb_metadata$mutant[x], comb_metadata$time.point[x], sep = "_")
  }
})

# Create cell-line specific SE data

## HCC3153
keep <- chip_oe_metadata$parental.cell.line == "HCC3153"
meta_hcc <- chip_oe_metadata[keep, ]
data_rpkm_norm_uf_hcc <- chip_oe_se_all_rpkm_norm_uf[, keep]
data_rpkm_uf_hcc <- chip_oe_se_all_rpkm[,keep]

# remove 0 variance rows
keep <- apply(data_rpkm_uf_hcc, 1, var) != 0
data_rpkm_hcc <- data_rpkm_uf_hcc[keep, ]

# filter low rpkm regions
keep <- apply(data_rpkm_hcc, 1, function(x) length(which(x > 1)) > 1)
data_rpkm_hcc <- data_rpkm_hcc[keep, ]

# normalize
data_rpkm_norm_hcc <- log2((data_rpkm_hcc*10) + 1)

## SUM185
keep <- chip_oe_metadata$parental.cell.line == "SUM185"
meta_sum <- chip_oe_metadata[keep, ]
data_rpkm_norm_uf_sum <- chip_oe_se_all_rpkm_norm_uf[, keep]
data_rpkm_uf_sum <- chip_oe_se_all_rpkm[,keep]

# remove 0 variance rows
keep <- apply(data_rpkm_uf_sum, 1, var) != 0
data_rpkm_sum <- data_rpkm_uf_sum[keep, ]

# filter low rpkm regions
keep <- apply(data_rpkm_sum, 1, function(x) length(which(x > 1)) > 1)
data_rpkm_sum <- data_rpkm_sum[keep, ]

# normalize
data_rpkm_norm_sum <- log2((data_rpkm_sum*10) + 1)

##########################################################################################
## Unsupervised analysis
##########################################################################################

pp_centered <- sweep(chip_oe_se_all_rpkm_norm, 1, rowMeans(chip_oe_se_all_rpkm_norm))

# subset to top 20% most variable features
chip_oe_se_ranks <- variable_features(pp_centered)
num_feats <- round(0.2 * length(chip_oe_se_ranks))
pp_var <- pp_centered[chip_oe_se_ranks[1:num_feats],]

# cluster
col_dist <- dist(t(pp_var), method = "euclidean")
col_hclust <- hclust(col_dist, method = "ward.D2")
row_dist <- dist(pp_var, method = "euclidean")
row_hclust <- hclust(row_dist, method = "ward.D2")

# form annotations
col_ann <- comb_metadata[,c("type", "mutant", "time.point", "condition")]
col_ann$type <- c("mesenchymal" = "mes", "basal" = "bas", "luminal" = "lum")[col_ann$type]
col_colors <- col_map
col_colors$type <- col_colors$type[c("bas", "lum")]

# cap
pp <- t(apply(pp_var, 1, function(x) sapply(x , function(y) min(max(-cap_val,y),cap_val))))

# plot
pdf("analysis/chip_seq_prrx1_oe/unsupervised/unsupervised_heatmap.pdf")
pheatmap(pp,
         annotation_col = col_ann,
         annotation_colors = col_colors,
         cluster_rows = row_hclust,
         cluster_cols = col_hclust,
         show_rownames = FALSE)
dev.off()

##########################################################################################
## Unsupervised analysis split by cell line
##########################################################################################

## HCC3153
pp_centered <- sweep(data_rpkm_norm_hcc, 1, rowMeans(data_rpkm_norm_hcc))

# subset to top 20% most variable features
chip_oe_se_ranks <- variable_features(pp_centered)
num_feats <- round(0.2 * length(chip_oe_se_ranks))
pp_var <- pp_centered[chip_oe_se_ranks[1:num_feats],]

# cluster
col_dist <- dist(t(pp_var), method = "euclidean")
col_hclust <- hclust(col_dist, method = "ward.D2")
row_dist <- dist(pp_var, method = "euclidean")
row_hclust <- hclust(row_dist, method = "ward.D2")

# form annotations
col_ann <- comb_metadata[,c("mutant", "time.point", "condition")]
col_colors <- col_map

# cap
pp <- t(apply(pp_var, 1, function(x) sapply(x , function(y) min(max(-cap_val,y),cap_val))))

# plot
pdf("analysis/chip_seq_prrx1_oe/unsupervised/unsupervised_heatmap_hcc3153.pdf")
pheatmap(pp,
         annotation_col = col_ann,
         annotation_colors = col_colors,
         cluster_rows = row_hclust,
         cluster_cols = col_hclust,
         show_rownames = FALSE)
dev.off()

## SUM185
pp_centered <- sweep(data_rpkm_norm_sum, 1, rowMeans(data_rpkm_norm_sum))

# subset to top 20% most variable features
chip_oe_se_ranks <- variable_features(pp_centered)
num_feats <- round(0.2 * length(chip_oe_se_ranks))
pp_var <- pp_centered[chip_oe_se_ranks[1:num_feats],]

# cluster
col_dist <- dist(t(pp_var), method = "euclidean")
col_hclust <- hclust(col_dist, method = "ward.D2")
row_dist <- dist(pp_var, method = "euclidean")
row_hclust <- hclust(row_dist, method = "ward.D2")

# form annotations
col_ann <- comb_metadata[,c("mutant", "time.point", "condition")]
col_colors <- col_map

# cap
pp <- t(apply(pp_var, 1, function(x) sapply(x , function(y) min(max(-cap_val,y),cap_val))))

# plot
pdf("analysis/chip_seq_prrx1_oe/unsupervised/unsupervised_heatmap_sum185.pdf")
pheatmap(pp,
         annotation_col = col_ann,
         annotation_colors = col_colors,
         cluster_rows = row_hclust,
         cluster_cols = col_hclust,
         show_rownames = FALSE)
dev.off()

##########################################################################################
## MOFA projection
##########################################################################################

mofa_data <- comb_se_data
rownames(mofa_data) <- paste0(rownames(mofa_data), "_H3K27ac")

data_list <- list(chip.se = as.matrix(mofa_data))
weight_list <- list(chip.se = chip_weights)
call_data <- prepare_mofa_calls(data_list, weight_list)
mofa_scores <- call_mofa_factor_scores(call_data)

pp <- cbind(mofa_scores, comb_metadata)

pdf("analysis/chip_seq_prrx1_oe/mofa/f2_f3.pdf")
ggplot(pp, aes(x = Factor2, y = Factor3, color = type, shape = group, fill = condition)) +
  geom_point() +
  xlab("Factor 2 Score") + ylab("Factor 3 Score") +
  scale_color_manual(values = line_colours$type) +
  scale_shape_manual(values = c("untreated" = 16, "sh_wt_short" = 21, "sh_dH3_short" = 22, "sh_wt_long" = 23, "sh_dH3_long" = 24)) +
  scale_fill_manual(values = c("no_dox" = "white", "dox" = "grey")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  geom_text_repel(aes(label = label), max.overlaps = Inf, show.legend = FALSE) +
  my_theme()
dev.off()

pdf("analysis/chip_seq_prrx1_oe/mofa/f2_f6.pdf")
ggplot(pp, aes(x = Factor2, y = Factor6, color = type, shape = group, fill = condition)) +
  geom_point() +
  xlab("Factor 2 Score") + ylab("Factor 6 Score") +
  scale_color_manual(values = line_colours$type) +
  scale_shape_manual(values = c("untreated" = 16, "sh_wt_short" = 21, "sh_dH3_short" = 22, "sh_wt_long" = 23, "sh_dH3_long" = 24)) +
  scale_fill_manual(values = c("no_dox" = "white", "dox" = "grey")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  geom_text_repel(aes(label = label), max.overlaps = Inf, show.legend = FALSE) +
  my_theme()
dev.off()

pdf("analysis/chip_seq_prrx1_oe/mofa/f3_f6.pdf")
ggplot(pp, aes(x = Factor3, y = Factor6, color = type, shape = group, fill = condition)) +
  geom_point() +
  xlab("Factor 3 Score") + ylab("Factor 6 Score") +
  scale_color_manual(values = line_colours$type) +
  scale_shape_manual(values = c("untreated" = 16, "sh_wt_short" = 21, "sh_dH3_short" = 22, "sh_wt_long" = 23, "sh_dH3_long" = 24)) +
  scale_fill_manual(values = c("no_dox" = "white", "dox" = "grey")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  geom_text_repel(aes(label = label), max.overlaps = Inf, show.legend = FALSE) +
  my_theme()
dev.off()

##########################################################################################
## Subtype SE TF's
##########################################################################################

# subset and rename rows
bas_se_ids <- rename_duplicates(chip_se_ann[bas_ses, "gene.symbol"])
mes_se_ids <- rename_duplicates(chip_se_ann[mes_ses, "gene.symbol"])
lum_se_ids <- rename_duplicates(chip_se_ann[lum_ses, "gene.symbol"])

pp_centered <- sweep(chip_oe_se_all_rpkm_norm_uf, 1, rowMeans(chip_oe_se_all_rpkm_norm_uf))
pp <- pp_centered[c(bas_ses, lum_ses, mes_ses),]
rownames(pp) <- c(bas_se_ids, lum_se_ids, mes_se_ids)

# cluster
col_dist <- dist(t(pp), method = "euclidean")
col_hclust <- hclust(col_dist, method = "ward.D2")
row_dist <- dist(pp, method = "euclidean")
row_hclust <- hclust(row_dist, method = "ward.D2")

# form annotations
col_ann <- comb_metadata[,c("type", "mutant", "time.point", "condition")]
col_ann$type <- c("mesenchymal" = "mes", "basal" = "bas", "luminal" = "lum")[col_ann$type]
col_colors <- col_map

se_annotation <- data.frame(row.names = rownames(pp),
                            type = c(rep("bas", length(bas_ses)),
                                     rep("lum", length(lum_ses)),
                                     rep("mes", length(mes_ses))))

# plot
pdf("analysis/chip_seq_prrx1_oe/subtype_tfs/subtype_tfs.pdf")
pheatmap(pp, 
         annotation_col = col_ann,
         annotation_row = se_annotation,
         annotation_colors = col_colors,
         cluster_rows = row_hclust,
         cluster_cols = col_hclust,
         show_rownames = TRUE)
dev.off()

##########################################################################################
## Subtype SE TF's by cell line
##########################################################################################

bas_ut_se_ids <- rename_duplicates(chip_se_ann[bas_ut_ses, "gene.symbol"])
mes_ut_se_ids <- rename_duplicates(chip_se_ann[mes_ut_ses, "gene.symbol"])
lum_ut_se_ids <- rename_duplicates(chip_se_ann[lum_ut_ses, "gene.symbol"])

# HCC3153

# subset and rename rows
pp_centered <- sweep(data_rpkm_norm_uf_hcc, 1, rowMeans(data_rpkm_norm_uf_hcc))
pp <- pp_centered[c(bas_ut_ses, lum_ut_ses, mes_ut_ses),]
rownames(pp) <- c(bas_ut_se_ids, lum_ut_se_ids, mes_ut_se_ids)

# cluster
col_dist <- dist(t(pp), method = "euclidean")
col_hclust <- hclust(col_dist, method = "ward.D2")
row_dist <- dist(pp, method = "euclidean")
row_hclust <- hclust(row_dist, method = "ward.D2")

# form annotations
col_ann <- comb_metadata[,c("mutant", "time.point", "condition")]
col_colors <- col_map

se_annotation <- data.frame(row.names = rownames(pp),
                            type = c(rep("bas", length(bas_ut_ses)),
                                     rep("lum", length(lum_ut_ses)),
                                     rep("mes", length(mes_ut_ses))))

# plot
pdf("analysis/chip_seq_prrx1_oe/subtype_tfs/subtype_ut_tfs_hcc3153.pdf",
    width = 6, height = 9)
pheatmap(pp, 
         annotation_col = col_ann,
         annotation_row = se_annotation,
         annotation_colors = col_colors,
         cluster_rows = row_hclust,
         cluster_cols = col_hclust,
         show_rownames = TRUE)
dev.off()

# SUM185

# subset and rename rows
pp_centered <- sweep(data_rpkm_norm_uf_sum, 1, rowMeans(data_rpkm_norm_uf_sum))
pp <- pp_centered[c(bas_ut_ses, lum_ut_ses, mes_ut_ses),]
rownames(pp) <- c(bas_ut_se_ids, lum_ut_se_ids, mes_ut_se_ids)

# cluster
col_dist <- dist(t(pp), method = "euclidean")
col_hclust <- hclust(col_dist, method = "ward.D2")
row_dist <- dist(pp, method = "euclidean")
row_hclust <- hclust(row_dist, method = "ward.D2")

# form annotations
col_ann <- comb_metadata[,c("mutant", "time.point", "condition")]
col_colors <- col_map

se_annotation <- data.frame(row.names = rownames(pp),
                            type = c(rep("bas", length(bas_ut_ses)),
                                     rep("lum", length(lum_ut_ses)),
                                     rep("mes", length(mes_ut_ses))))

# plot
pdf("analysis/chip_seq_prrx1_oe/subtype_tfs/subtype_ut_tfs_sum185.pdf",
    width = 6, height = 9)
pheatmap(pp, 
         annotation_col = col_ann,
         annotation_row = se_annotation,
         annotation_colors = col_colors,
         cluster_rows = row_hclust,
         cluster_cols = col_hclust,
         show_rownames = TRUE)
dev.off()
