library(Seurat)
library(foreach)
library(doMC)
registerDoMC(26)

diff_test <- function(sc){
  clust_ids <- sort(as.numeric(as.character(unique(sc@meta.data$seurat_clusters))))
  diff_list <- list()
  for(i in 1:length(clust_ids)){
    c_id <- clust_ids[i]
    message("Cluster ", c_id, " of ", max(clust_ids))
    diff_list[[i]] <- FindMarkers(sc, ident.1 = c_id, assay = "RNA", slot = "data")
  }
  return(list(clust_ids, diff_list))
}

# Load data ------------------------------------------------------------------------------------------

setwd("../")

load(file = "R_Data/SC_List_Ind_Norm_Clust.RData")
load(file = "R_Data/SC_Metadata.RData")

# Differential expression ----------------------------------------------------------------------------

res <- foreach(i = 1:length(sc_list)) %dopar% {
  diff_test(sc_list[[i]])
}
names(res) <- scrna_sample_table$sample3
save(res, file = "R_Data/SC_Clust_DE_Genes.RData")

for(i in 1:length(res)){
  sample <- names(res)[i]
  clust_ids <- res[[i]][[1]]
  df_list <- res[[i]][[2]]
  df_list <- lapply(1:length(df_list), function(x){
    df <- cbind(gene = rownames(df_list[[x]]), df_list[[x]], clust = clust_ids[x])
    rownames(df) <- NULL
    df[which(df$p_val_adj < 0.05),]
  })
  out <- do.call(rbind, df_list)
  fn <- paste0("analysis/scrna_seq/cell_lines/de_genes/", sample, "_de_genes.csv")
  write.table(out, file = fn, sep = ",", row.names = FALSE)
}

## save processed output
de_genes <- list()
for(i in 1:length(res)){
  sample <- names(res)[i]
  clust_ids <- res[[i]][[1]]
  df_list <- res[[i]][[2]]
  df_list <- lapply(1:length(df_list), function(x){
    df <- cbind(gene = rownames(df_list[[x]]), df_list[[x]], clust = clust_ids[x])
    rownames(df) <- NULL
    df[which(df$p_val_adj < 0.05),]
  })
  de_genes[[i]] <- do.call(rbind, df_list)
}
names(de_genes) <- names(res)
save(de_genes, file = "R_Data/SC_DE_Genes.RData")
