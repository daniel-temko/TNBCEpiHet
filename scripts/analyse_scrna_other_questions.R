setwd("../")

library(Seurat)

load(file = "R_Data/SC_BC_DR.RData")
load(file = "R_Data/Heatmap_Metadata.RData")

sc$parental.tnbc.type <- line_annotation[sc$orig.ident, "type"]

pdf("analysis/scrna_seq/other_analyses/igfbp7.pdf")
p <- FeaturePlot(sc, "IGFBP7")
print(p)
dev.off()

pdf("analysis/scrna_seq/other_analyses/cd44.pdf")
p <- FeaturePlot(sc, "CD44")
print(p)
dev.off()

Idents(sc) <- "parental.tnbc.type"
res <- FindMarkers(sc, ident.1 = "mesenchymal")
