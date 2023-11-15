library(pheatmap)
source("../utilities/utils.R")

setwd("../")

up_tfs <- function(df){
  df[which((df$log2FoldChange > 0) & (df$dbd.tf)),]
}

up_ut_tfs <- function(df, tf_col){
  df[which((df$log2FoldChange > 0) & (df$ut.tf)),]
}

load("R_Data/RNA_Data.RData")
load("R_Data/Heatmap_Metadata.RData")

bas_de <- read_csv("analysis/rna_seq/differential_genes/basal_de_lfc0_genes.csv")
mes_de <- read_csv("analysis/rna_seq/differential_genes/mesenchymal_de_lfc0_genes.csv")
lum_de <- read_csv("analysis/rna_seq/differential_genes/luminal_de_lfc0_genes.csv")

bas_dse <- read_csv("analysis/chip_seq/differential_regions/superenhancers_basal_diff_lfc0.csv")
mes_dse <- read_csv("analysis/chip_seq/differential_regions/superenhancers_mesenchymal_diff_lfc0.csv")
lum_dse <- read_csv("analysis/chip_seq/differential_regions/superenhancers_luminal_diff_lfc0.csv")

bas_tfs <- intersect(up_tfs(bas_de)$gene.symbol, up_tfs(bas_dse)$closest.gene)
mes_tfs <- intersect(up_tfs(mes_de)$gene.symbol, up_tfs(mes_dse)$closest.gene)
lum_tfs <- intersect(up_tfs(lum_de)$gene.symbol, up_tfs(lum_dse)$closest.gene)

# Get UT TFs
bas_ut_tfs <- intersect(up_ut_tfs(bas_de)$gene.symbol, up_ut_tfs(bas_dse)$closest.gene)
mes_ut_tfs <- intersect(up_ut_tfs(mes_de)$gene.symbol, up_ut_tfs(mes_dse)$closest.gene)
lum_ut_tfs <- intersect(up_ut_tfs(lum_de)$gene.symbol, up_ut_tfs(lum_dse)$closest.gene)

subtype_tf_df <- data.frame(tf = c(bas_ut_tfs, mes_ut_tfs, lum_ut_tfs),
                            type = c(rep("basal", length(bas_ut_tfs)),
                                     rep("mesenchymal", length(mes_ut_tfs)),
                                     rep("luminal", length(lum_ut_tfs))))

# Output subtype TF's to file
write.table(subtype_tf_df,
            file = "analysis/integration/subtype_regulation/subtype_tfs.csv",
            sep = ",",
            col.names = NA)

data <- rna_rpkm_norm_uf[c(bas_tfs, mes_tfs, lum_tfs), ]
col_types <- line_annotation[colnames(data), "type"]
data <- data[,order(col_types)]
data <- sweep(data, 1, rowMeans(data))

pdf("analysis/integration/subtype_regulation/de_dse_tfs_lfc.pdf")
  p <- pheatmap(data, 
           annotation_col = line_annotation, 
           annotation_colors = line_colours, 
           cluster_rows = FALSE,
           cluster_cols = FALSE)
  print(p)
dev.off()

data_std <- t(apply(data, 1, function(x) x / sd(x)))

pdf("analysis/integration/subtype_regulation/de_dse_tfs_row_norm.pdf")
  p <- pheatmap(data_std, 
                annotation_col = line_annotation, 
                annotation_colors = line_colours, 
                cluster_rows = FALSE,
                cluster_cols = FALSE)
  print(p)
dev.off()

# row ordered version
data <- rna_rpkm_norm_uf[c(sort(bas_tfs), sort(mes_tfs), sort(lum_tfs)), ]
col_types <- line_annotation[colnames(data), "type"]
data <- data[,order(col_types)]
data <- sweep(data, 1, rowMeans(data))

pdf("analysis/integration/subtype_regulation/de_dse_tfs_lfc_row_ordered.pdf")
  p <- pheatmap(data, 
           annotation_col = line_annotation, 
           annotation_colors = line_colours, 
           cluster_rows = FALSE,
           cluster_cols = FALSE)
  print(p)
dev.off()

# clustered version
data <- rna_rpkm_norm_uf[c(bas_tfs, mes_tfs, lum_tfs), ]
data <- sweep(data, 1, rowMeans(data))

col_dist <- dist(t(data), method = "euclidean")
col_hclust <- hclust(col_dist, method = "ward.D2")

row_dist <- dist(data, method = "euclidean")
row_hclust <- hclust(row_dist, method = "ward.D2")

tf_annotation <- data.frame(row.names = rownames(data), 
                            type = c(rep("basal", length(bas_tfs)),
                                     rep("mesenchymal", length(mes_tfs)),
                                     rep("luminal", length(lum_tfs))))

pdf("analysis/integration/subtype_regulation/de_dse_tfs_lfc_clust.pdf")
  p <- pheatmap(data, 
           annotation_col = line_annotation, 
           annotation_row = tf_annotation,
           annotation_colors = line_colours, 
           cluster_rows = row_hclust,
           cluster_cols = col_hclust)
dev.off()

# UT clustered version
data_ut <- rna_rpkm_norm_uf[c(bas_ut_tfs, mes_ut_tfs, lum_ut_tfs), ]
data_ut <- sweep(data_ut, 1, rowMeans(data_ut))

col_dist_ut <- dist(t(data_ut), method = "euclidean")
col_hclust_ut <- hclust(col_dist_ut, method = "ward.D2")

row_dist_ut <- dist(data_ut, method = "euclidean")
row_hclust_ut <- hclust(row_dist_ut, method = "ward.D2")

tf_annotation_ut <- data.frame(row.names = rownames(data_ut), 
                               type = c(rep("basal", length(bas_ut_tfs)),
                                        rep("mesenchymal", length(mes_ut_tfs)),
                                        rep("luminal", length(lum_ut_tfs))))

pdf("analysis/integration/subtype_regulation/de_dse_tfs_lfc_clust_ut.pdf")
p <- pheatmap(data_ut, 
              annotation_col = line_annotation, 
              annotation_row = tf_annotation_ut,
              annotation_colors = line_colours, 
              cluster_rows = row_hclust_ut,
              cluster_cols = col_hclust_ut)
dev.off()
