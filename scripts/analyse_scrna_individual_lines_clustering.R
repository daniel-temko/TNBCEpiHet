library(Seurat)
library(SingleCellExperiment)

# Load data ------------------------------------------------------------------------------------------

setwd("../")

load(file = "R_Data/SC_List_Ind_Norm.RData")
load(file = "R_Data/SC_Metadata.RData")

# Analysis ------------------------------------------------------------------------------------------

# Clustering
for(i in 1:length(sc_list)){
  sc_list[[i]] <- FindVariableFeatures(sc_list[[i]])
  sc_list[[i]] <- RunPCA(sc_list[[i]])
  sc_list[[i]] <- RunUMAP(sc_list[[i]], dims = 1:20)
  
  sc_list[[i]] <- FindNeighbors(sc_list[[i]], dims = 1:20)
  sc_list[[i]] <- FindClusters(sc_list[[i]], resolution = 0.5)
}
save(sc_list, file = "R_Data/SC_List_Ind_Norm_Clust.RData")
