library(DESeq2)

setwd("../")

load("R_Data/PDX_RNA_Data.RData")

sample.sheet <- data.frame(sample = colnames(pdx_rna_counts))
dds <- DESeqDataSetFromMatrix(countData = pdx_rna_counts,
                              colData = sample.sheet,
                              design = ~ 1)
dds <- estimateSizeFactors(dds)
pdx_rna_counts_norm <- counts(dds, normalized = TRUE)

write.table(pdx_rna_counts, file = "analysis/pdx_rna_seq/summary/pdx_rna_counts.csv", sep = ",", col.names = NA)
write.table(pdx_rna_counts_norm, file = "analysis/pdx_rna_seq/summary/pdx_rna_counts_norm.csv", sep = ",", col.names = NA)
write.table(pdx_rna_rpkm_norm, file = "analysis/pdx_rna_seq/summary/pdx_rna_fpkm_norm.csv", sep = ",", col.names = NA)
write.table(pdx_rna_rpkm_norm_uf, file = "analysis/pdx_rna_seq/summary/pdx_rna_fpkm_norm_unfiltered.csv", sep = ",", col.names = NA)
