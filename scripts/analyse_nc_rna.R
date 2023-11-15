source("../utilities/preprocessing.R")
source("../utilities/utils.R")
library(pheatmap)
library(ggplot2)
library(ggrepel)

setwd("../")

## Local functions ------------------------------------------------------------------------------------------------------------------------------------

up_tfs <- function(df){
  df[which((df$log2FoldChange > 0) & (df$dbd.tf)),]
}

## Load data ------------------------------------------------------------------------------------------------------------------------------------

load("R_Data/RNA_NC_Combined_Data.RData")
load("R_Data/RNA_Data.RData")

bas_de <- read_csv("analysis/rna_seq/differential_genes/basal_de_lfc0_genes.csv")
mes_de <- read_csv("analysis/rna_seq/differential_genes/mesenchymal_de_lfc0_genes.csv")
lum_de <- read_csv("analysis/rna_seq/differential_genes/luminal_de_lfc0_genes.csv")

bas_dse <- read_csv("analysis/chip_seq/differential_regions/superenhancers_basal_diff_lfc0.csv")
mes_dse <- read_csv("analysis/chip_seq/differential_regions/superenhancers_mesenchymal_diff_lfc0.csv")
lum_dse <- read_csv("analysis/chip_seq/differential_regions/superenhancers_luminal_diff_lfc0.csv")

bas_tfs <- intersect(up_tfs(bas_de)$gene.symbol, up_tfs(bas_dse)$closest.gene)
mes_tfs <- intersect(up_tfs(mes_de)$gene.symbol, up_tfs(mes_dse)$closest.gene)
lum_tfs <- intersect(up_tfs(lum_de)$gene.symbol, up_tfs(lum_dse)$closest.gene)
reg_tfs <- c(bas_tfs, mes_tfs, lum_tfs)
reg_tfs_ord <- c(sort(bas_tfs), sort(mes_tfs), sort(lum_tfs))

colour_code <- read.table("metadata/colour_code_coclustering.csv", sep = ",", header = TRUE, comment.char = "")
colour_map <- as.character(colour_code$colour)
names(colour_map) <- colour_code$type

prop_features <- 0.2

rna_comb_metadata$type <- with(rna_comb_metadata, sapply(1:nrow(rna_comb_metadata), function(x) {
  if(cancer.type[x] == "breast") {
    if(breast.cancer.subtype[x] == "tnbc") {
      paste0("breast.", tnbc.subtype[x])  
    } else {
      "breast.other"  
    }
  } else {
    cancer.type[x]
  }
}))

## PRRX1 expression ------------------------------------------------------------------------------------------------------------------------------------

nn <- c("Hs578", "MDAMB157", "MDAMB436", "UACC3199", "HDQP1", "R_TTC642", "R_GIMEN", "R_LS")
ee <- rna_comb_counts_norm["PRRX1", nn, drop = FALSE]
ll <- nn
ll[c(1, 6, 7, 8)] <- c("Hs578T", "TTC642", "GIMEN", "LS")
colnames(ee) <- ll
pdf("analysis/rna_seq_nc/prrx1_ln_exp.pdf", 6, 5)
p <- pheatmap(ee, cluster_rows = FALSE, cluster_cols = FALSE, cellheight = 40)
print(p)
dev.off()

## Counts ------------------------------------------------------------------------------------------------------------------------------------

n_features <- round(prop_features * nrow(rna_comb_counts_norm))
var_feats <- variable_features(rna_comb_counts_norm)

# Heatmap
data <- rna_comb_counts_norm[var_feats[1:n_features],]

data <- sweep(data, 1, rowMeans(data))

col_dist <- dist(t(data), method = "euclidean")
col_hclust <- hclust(col_dist, method = "ward.D2")
row_dist <- dist(data, method = "euclidean")
row_hclust <- hclust(row_dist, method = "ward.D2")

data <- t(apply(data, 1, function(x) sapply(x , function(y) min(max(-6,y),6))))

annotation_col <- data.frame(row.names = rna_comb_metadata$sample,
                             type = rna_comb_metadata$type,
                             source = rna_comb_metadata$source)
annotation_colors <- list(type = colour_map)

pdf("analysis/rna_seq_nc/counts_norm_top_20pc_combined_heatmap.pdf", width = 11, height = 9)
p <- pheatmap(data,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              cluster_rows = row_hclust,
              cluster_cols = col_hclust,
              show_rownames = FALSE)
print(p)
dev.off()

# PCA filtering
data <- rna_comb_counts_norm[var_feats[1:n_features],]

data_pca <- prcomp(t(as.matrix(data)))

data <- as.data.frame(data_pca$x)[,c("PC1", "PC2")]
data$cell.line <- as.character(rna_comb_metadata$cell.line)
data$type <- as.character(rna_comb_metadata$type)
data$label <- rownames(data)
data$label <- sapply(1:nrow(data), function(x){
  if((data$type[x] %in% c("breast.mesenchymal", "rhabdoid")) | (data$cell.line[x] == "GIMEN")) {
    data$label[x]
  } else {
    ""
  }
})
data$source <- as.character(rna_comb_metadata$source)
pdf(file = paste0("analysis/rna_seq_nc/counts_norm_top_20pc_pca.pdf"))
p <- ggplot(data = data, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = type, shape = source)) +
  scale_color_manual(values = c(colour_map)) +
  geom_text_repel(aes(label = label), max.overlaps = Inf) +
  theme(legend.position = "top",
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16)) +
  guides(color = guide_legend(ncol = 2, title.position = "top"), 
         shape = guide_legend(ncol = 1, title.position = "top"))
print(p)
dev.off()

# SE heatmap - Note: unfiltered data is not needed here since all TF's are present in filtered
# data
stopifnot(all(reg_tfs %in% rownames(rna_comb_counts_norm)))

data <- rna_comb_counts_norm[reg_tfs, ]
data <- sweep(data, 1, rowMeans(data))

annotation_col <- data.frame(row.names = rna_comb_metadata$sample,
                             type = rna_comb_metadata$type,
                             source = rna_comb_metadata$source)
annotation_colors <- list(type = colour_map)

data <- data[,order(annotation_col$type)]

pdf("analysis/rna_seq_nc/count_norm_se_tfs_combined_heatmap.pdf", width = 12, height = 8)
par(mar = c(10,4,4,4))
p <- pheatmap(data,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              show_rownames = TRUE)
print(p)
dev.off()

data_std <- t(apply(data, 1, function(x) x / sd(x)))

pdf("analysis/rna_seq_nc/count_norm_se_tfs_combined_heatmap_row_norm.pdf", width = 12, height = 8)
par(mar = c(10,4,4,4))
p <- pheatmap(data_std,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              show_rownames = TRUE)
print(p)
dev.off()

# row ordered
data <- rna_comb_counts_norm[reg_tfs_ord, ]
data <- sweep(data, 1, rowMeans(data))

annotation_col <- data.frame(row.names = rna_comb_metadata$sample,
                             type = rna_comb_metadata$type,
                             source = rna_comb_metadata$source)
annotation_colors <- list(type = colour_map)

data <- data[,order(annotation_col$type)]

pdf("analysis/rna_seq_nc/count_norm_se_tfs_combined_heatmap_row_ordered.pdf", width = 12, height = 8)
par(mar = c(10,4,4,4))
p <- pheatmap(data,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              show_rownames = TRUE)
print(p)
dev.off()

# Clustered version
data <- rna_comb_counts_norm[reg_tfs, ]
data <- sweep(data, 1, rowMeans(data))

col_dist <- dist(t(data), method = "euclidean")
col_hclust <- hclust(col_dist, method = "ward.D2")

row_dist <- dist(data, method = "euclidean")
row_hclust <- hclust(row_dist, method = "ward.D2")

annotation_col <- data.frame(row.names = rna_comb_metadata$sample,
                             type = rna_comb_metadata$type,
                             source = rna_comb_metadata$source)
annotation_row <- data.frame(row.names = rownames(data),
                             type = c(rep("breast.basal", length(bas_tfs)),
                                      rep("breast.mesenchymal", length(mes_tfs)),
                                      rep("breast.luminal", length(lum_tfs))))
annotation_colors <- list(type = colour_map)

pdf("analysis/rna_seq_nc/count_norm_se_tfs_combined_heatmap_clust.pdf", width = 12, height = 8)
par(mar = c(10,4,4,4))
p <- pheatmap(data,
              annotation_col = annotation_col,
              annotation_row = annotation_row,
              annotation_colors = annotation_colors,
              cluster_rows = row_hclust,
              cluster_cols = col_hclust)
print(p)
dev.off()

## RPKM ------------------------------------------------------------------------------------------------------------------------------------

# HVG selection
n_features <- round(prop_features * nrow(rna_comb_rpkm_norm))
var_feats <- variable_features(rna_comb_rpkm_norm)

# Plot data
data <- rna_comb_rpkm_norm[var_feats[1:n_features],]

data <- sweep(data, 1, rowMeans(data))

col_dist <- dist(t(data), method = "euclidean")
col_hclust <- hclust(col_dist, method = "ward.D2")
row_dist <- dist(data, method = "euclidean")
row_hclust <- hclust(row_dist, method = "ward.D2")

data <- t(apply(data, 1, function(x) sapply(x , function(y) min(max(-6,y),6))))

annotation_col <- data.frame(row.names = rna_comb_metadata$sample,
                             type = rna_comb_metadata$type,
                             source = rna_comb_metadata$source)
annotation_colors <- list(type = colour_map)

pdf("analysis/rna_seq_nc/rpkm_norm_top_20pc_combined_heatmap.pdf", width = 11.69, height = 8.27)
p <- pheatmap(data,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              cluster_rows = row_hclust,
              cluster_cols = col_hclust,
              show_rownames = FALSE)
print(p)
dev.off()

# PCA filtering
data <- rna_comb_rpkm_norm[var_feats[1:n_features],]

data_pca <- prcomp(t(as.matrix(data)))

type <- as.character(rna_comb_metadata$type)
source <- as.character(rna_comb_metadata$source)
data <- as.data.frame(data_pca$x)
pdf(file = paste0("analysis/rna_seq_nc/rpkm_norm_top_20pc_pca.pdf"), width = 11.69, height = 8.27)
p <- ggplot(data = data, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = type, shape = source)) +
  geom_text_repel(label = rna_comb_metadata$cell.line, max.overlaps = Inf) +
  scale_color_manual(values = c(colour_map))
print(p)
dev.off()
