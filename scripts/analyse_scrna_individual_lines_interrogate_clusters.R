source("../utilities/utils.R")
source("../utilities/sampling.R")
source("../utilities/signatures.R")
library(Seurat)
library(ggplot2)
library(reshape2)

# Load data ------------------------------------------------------------------------------------------

setwd("../")

load(file = "R_Data/SC_List_Ind_Norm_Clust.RData")
load(file = "R_Data/SC_Metadata.RData")
load(file = "R_Data/Gene_Sigs.RData")
load(file = "R_Data/TF_Data.RData")

# Plot clustering ------------------------------------------------------------------------------------------------

for(i in 1:length(sc_list)){
  sample <- scrna_sample_table$sample3[i]
  f_name <- paste0("analysis/scrna_seq/cell_lines/clustering/", sample, "_clusters.pdf")
  pdf(f_name)
  p <- DimPlot(sc_list[[i]], label = TRUE) + NoLegend()  
  print(p)
  dev.off()
}

# Analyze cluster differences ------------------------------------------------------------------------------------

# TNBC type signature scores
tnbc_sigs <- list()
for(i in 1:length(sc_list)){
  message("\n", i)
  exp_data <- as.matrix(sc_list[[i]][["RNA"]]@data)
  exp_data <- sweep(exp_data, 1, rowMeans(exp_data))
  ids <- BalancedSample(labels = sc_list[[i]]$seurat_clusters, size = 100, seed = 20210113)
  exp_data <- exp_data[,ids]
  clust_vals <- sc_list[[i]]$seurat_clusters[ids]
  
  tnbc_gsl_sub <- uniform_subset_gsl(tnbc_gsl, rownames(exp_data))
  tnbc_sig_scores <- score_gsl(exp_matrix = exp_data, gsl = tnbc_gsl_sub, Nthreads = 32)
  tnbc_sigs[[i]] <- cbind(tnbc_sig_scores, cluster = clust_vals)
}
save(tnbc_sigs, file = "R_Data/SC_Sig_Scores.RData")

# Violin plots
for(i in 1:length(sc_list)){
  sample <- scrna_sample_table$sample3[i]
  data <- melt(tnbc_sigs[[i]])
  
  pdf(paste0("analysis/scrna_seq/cell_lines/signatures/", sample, "_tnbc_sig.pdf"))
  p <- ggplot(data, aes(x = cluster, y = value, group = cluster, fill = cluster)) + 
    geom_violin() +
    facet_grid(. ~ variable)
  print(p)
  dev.off()
}

# HDQP1 sub-cluster analysis
id <- match("HDQP1", scrna_sample_table$cell.line)
diff_genes <- FindMarkers(sc_list[[id]], ident.1 = 6)
diff_genes$dbd.tf <- rownames(diff_genes) %in% tf_df$Gene.name
diff_genes$ut.tf <- rownames(diff_genes) %in% tf_ut$HGNC.symbol

write.table(diff_genes, "analysis/scrna_seq/cell_lines/hdqp1/gene_lists/clust_6_genes.csv", 
            sep = ",",
            col.names = NA)

diff_tfs <- diff_genes[diff_genes$dbd.tf & (diff_genes$avg_logFC > 0),]

pdf("analysis/scrna_seq/cell_lines/hdqp1/plots/diff_tfs.pdf")
DotPlot(sc_list[[id]], features = rownames(diff_tfs)) + 
  RotatedAxis()
dev.off()

pdf("analysis/scrna_seq/cell_lines/hdqp1/plots/top_diff_tfs_umap.pdf")
FeaturePlot(sc_list[[id]], c("NR2F1", "PRRX1", "RUNX3", "ETV1"))
dev.off()

pdf("analysis/scrna_seq/cell_lines/hdqp1/plots/top6_diff_tfs_umap.pdf", height = 11, width = 9)
FeaturePlot(sc_list[[id]], rownames(diff_tfs)[1:6])
dev.off()

diff_ut_tfs <- diff_genes[diff_genes$ut.tf & (diff_genes$avg_logFC > 0),]

pdf("analysis/scrna_seq/cell_lines/hdqp1/plots/diff_ut_tfs.pdf")
DotPlot(sc_list[[id]], features = rownames(diff_ut_tfs)) + 
  RotatedAxis()
dev.off()

pdf("analysis/scrna_seq/cell_lines/hdqp1/plots/top6_diff_ut_tfs_umap.pdf", height = 11, width = 9)
FeaturePlot(sc_list[[id]], rownames(diff_ut_tfs)[1:6])
dev.off()
