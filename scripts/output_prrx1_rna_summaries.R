library(DESeq2)

setwd("../")

load("R_Data/RNA_PRRX1_Data.RData")
load("R_Data/RNA_PRRX1_LT_Data.RData")

sample.sheet <- data.frame(sample = factor(colnames(rna_prrx1_counts)))
prrx1.dds <- DESeqDataSetFromMatrix(countData = rna_prrx1_counts, 
                                    colData = sample.sheet, 
                                    design = ~ sample)
prrx1.dds <- estimateSizeFactors(prrx1.dds)
rna_prrx1_counts_norm <- counts(prrx1.dds, normalized = TRUE)

write.table(rna_prrx1_counts, file  = "analysis/rna_seq_prrx1_kd/summary/rna_prrx1_counts.csv", sep = ",", col.names = NA)
write.table(rna_prrx1_counts_norm, file = "analysis/rna_seq_prrx1_kd/summary/rna_prrx1_counts_norm.csv", sep = ",", col.names = NA)
write.table(rna_prrx1_rpkm_norm, file = "analysis/rna_seq_prrx1_kd/summary/rna_prrx1_fpkm.csv", sep = ",", col.names = NA)
write.table(rna_prrx1_metadata, file = "analysis/rna_seq_prrx1_kd/summary/rna_prrx1_metadata.csv", sep = ",", col.names = NA)

sample.sheet <- data.frame(sample = factor(colnames(rna_lt_counts)))
prrx1.lt.dds <- DESeqDataSetFromMatrix(countData = rna_lt_counts,
                                       colData = sample.sheet,
                                       design = ~ sample)
prrx1.lt.dds <- estimateSizeFactors(prrx1.lt.dds)
rna_lt_counts_norm <- counts(prrx1.lt.dds, normalized = TRUE)

write.table(rna_lt_counts, file = "analysis/rna_seq_prrx1_lt/summary/rna_prrx1_lt_counts.csv", sep = ",", col.names = NA)
write.table(rna_lt_counts_norm, file = "analysis/rna_seq_prrx1_lt/summary/rna_prrx1_lt_counts_norm.csv", sep = ",", col.names = NA)
write.table(rna_lt_rpkm_norm, file = "analysis/rna_seq_prrx1_lt/summary/rna_prrx1_lt_fpkm.csv", sep = ",", col.names = NA)
write.table(rna_lt_metadata, file = "analysis/rna_seq_prrx1_lt/summary/rna_prrx1_lt_metadata.csv", sep = ",", col.names = NA)
