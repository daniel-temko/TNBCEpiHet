
setwd("../")

load("R_Data/Heatmap_Metadata.RData")
load("R_Data/RNA_Data.RData")

rna_metadata$tnbc.subtype <- factor(rna_metadata$tnbc.subtype, levels = c("basal", "luminal", "mesenchymal"))

pdf("analysis/rna_seq/other_analyses/igfbp7.pdf")
boxplot(as.numeric(rna_rpkm_norm["IGFBP7", ]) ~ rna_metadata$tnbc.subtype,
        col = line_colours$type[c("basal", "luminal", "mesenchymal")])
dev.off()