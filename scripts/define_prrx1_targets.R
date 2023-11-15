# Script to define PRRX1 transcriptional targets in mes lines using ChIP and know-down RNA-seq
# data
library(DESeq2)
library(stringi)
setwd("../")

# Load data ------------------------------------------------------------------------------------------------

load("R_Data/PRRX1_ChIP_Data.RData")
load("R_Data/RNA_PRRX1_DE_Genes.RData")

hs578_all <- as.data.frame(hs578_de_lfc0_genes)
tt642_all <- as.data.frame(tt642_de_lfc0_genes)

# ChIP-seq Hs578 only
hs578_targets <- rownames(prrx1_chip_cl_targets)[prrx1_chip_cl_targets$Hs578]
tt642_targets <- rownames(prrx1_chip_cl_targets)[prrx1_chip_cl_targets$TT642]

# ChIP-seq All mes lines
mes_lines <- prrx1_chip_metadata[which(prrx1_chip_metadata$cell.line.type == "mesenchymal"),]
ids <- apply(prrx1_chip_cl_targets[,unique(mes_lines$cell.line)], 1, sum) == 3
mes_targets <- rownames(prrx1_chip_cl_targets)[ids]

prrx1_gml <- list(hs578_targets = hs578_targets,
                  tt642_targets = tt642_targets,
                  mes_targets = mes_targets)

# Hs578 DE genes (Note reversal of up/dn due to KD experiment)
hs578_dn <- rownames(hs578_all)[which(hs578_all$log2FoldChange > 0)]
hs578_up <- rownames(hs578_all)[which(hs578_all$log2FoldChange < 0)]

# TT642 DE genes (Note reversal of up/dn due to KD experiment)
tt642_dn <- rownames(tt642_all)[which(tt642_all$log2FoldChange > 0)]
tt642_up <- rownames(tt642_all)[which(tt642_all$log2FoldChange < 0)]

prrx1_gsl <- list(hs578_rna_targets = list(up = intersect(hs578_up, prrx1_gml$hs578),
                                           dn = intersect(hs578_dn, prrx1_gml$hs578)),
                  tt642_rna_targets = list(up = intersect(tt642_up, prrx1_gml$tt642_targets),
                                           dn = intersect(tt642_dn, prrx1_gml$tt642_targets)),
                  mes_rna_targets = list(up = intersect(hs578_up, prrx1_gml$mes),
                                         dn = intersect(hs578_dn, prrx1_gml$mes)))

x <- stri_list2matrix(prrx1_gsl$mes_rna_targets)
colnames(x) <- c("up", "dn")
write.table(x, "analysis/prrx1_chip_seq/signatures/mes_rna_targets.csv", 
            sep = ",",
            row.names = FALSE)

save(prrx1_gml,
     prrx1_gsl,
     file = "R_Data/PRRX1_Targets.RData")