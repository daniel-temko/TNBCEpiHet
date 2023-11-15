library(DESeq2)

setwd("../")

load("R_Data/RNA_Data.RData")
load("R_Data/RNA_Metadata.RData")

sample.sheet <- data.frame(sample = colnames(rna_counts))
dds <- DESeqDataSetFromMatrix(countData =  rna_counts, 
                              colData = sample.sheet, 
                              design = ~ 1)
dds <- estimateSizeFactors(dds)
rna_counts_norm <- counts(dds, normalized = TRUE)

write.table(rna_counts, file = "analysis/rna_seq/summary/rna_counts.csv", sep = ",", col.names = NA)
write.table(rna_counts_norm, file = "analysis/rna_seq/summary/rna_counts_norm.csv", sep = ",", col.names = NA)
write.table(rna_rpkm_norm, file = "analysis/rna_seq/summary/rna_fpkm_norm.csv", sep = ",", col.names = NA)
write.table(rna_rpkm_norm_uf, file = "analysis/rna_seq/summary/rna_fpkm_norm_unfiltered.csv", sep = ",", col.names = NA)
