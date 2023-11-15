library(DESeq2)
library(ggplot2)
library(ggrepel)

setwd("../")

# Load data -------------------------------------------------------------------------------------------

load("R_Data/RNA_PRRX1_Data.RData")

# DESeq2 -------------------------------------------------------------------------------------------

col_data <- rna_prrx1_metadata

col_data$treat <- ! col_data$shrna == "no_sh"
col_data$induced <- col_data$treat & col_data$condition == "DOX"
col_data$group <- paste0(col_data$parental.cell.line, "_", col_data$induced)
col_data$group <- factor(col_data$group)

# Hs578 -------------------------------------------------------------------------------------------

col_data_sub <- droplevels(subset(col_data, parental.cell.line == "Hs578"))
counts_sub <- rna_prrx1_counts[,rownames(col_data_sub)]

dds <- DESeqDataSetFromMatrix(countData = counts_sub,
                              colData = col_data_sub,
                              design = ~ group)

keep <- rowSums(counts(dds)) > 0
dds <- dds[keep, ]

hs578_dds <- DESeq(dds)

hs578_res <- results(hs578_dds, contrast = c("group", "Hs578_TRUE", "Hs578_FALSE"))

hs578_sub <- subset(hs578_res, padj < 0.05)
hs578_sub <- hs578_sub[order(hs578_sub$padj), ]

write.table(hs578_sub, 
            file = "analysis/rna_seq_prrx1_kd/differential_genes/hs578_de_lfc0_genes.csv", 
            sep = ",",
            col.names = NA)
hs578_de_lfc0_genes <- hs578_sub

hs578_res <- results(hs578_dds, contrast = c("group", "Hs578_FALSE", "Hs578_TRUE"))

hs578_beta <- data.frame(gene = rownames(hs578_res), lfc = hs578_res$log2FoldChange, p.value = hs578_res$pvalue)
keep <- !is.na(hs578_beta$p.value)
hs578_beta <- hs578_beta[keep, ]
write.table(hs578_beta, 
            file = "analysis/prrx1_beta/hs578/hs578_beta", 
            sep = "\t", 
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

# TT642 -------------------------------------------------------------------------------------------

col_data_sub <- droplevels(subset(col_data, parental.cell.line == "TT642"))
counts_sub <- rna_prrx1_counts[,rownames(col_data_sub)]

dds <- DESeqDataSetFromMatrix(countData = counts_sub,
                              colData = col_data_sub,
                              design = ~ group)

keep <- rowSums(counts(dds)) > 0
dds <- dds[keep, ]

tt642_dds <- DESeq(dds)

tt642_res <- results(tt642_dds, contrast = c("group", "TT642_TRUE", "TT642_FALSE"))

tt642_sub <- subset(tt642_res, padj < 0.05)
tt642_sub <- tt642_sub[order(tt642_sub$padj), ]

write.table(tt642_sub, 
            file = "analysis/rna_seq_prrx1_kd/differential_genes/tt642_de_lfc0_genes.csv", 
            sep = ",",
            col.names = NA)
tt642_de_lfc0_genes <- tt642_sub

tt642_res <- results(tt642_dds, contrast = c("group", "TT642_FALSE", "TT642_TRUE"))

tt642_beta <- data.frame(gene = rownames(tt642_res), lfc = tt642_res$log2FoldChange, p.value = tt642_res$pvalue)
keep <- !is.na(tt642_beta$p.value)
tt642_beta <- tt642_beta[keep, ]
write.table(tt642_beta, 
            file = "analysis/prrx1_beta/tt642/tt642_beta", 
            sep = "\t", 
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

# Save results --------------------------------------------------------------------------------------------------------

save(hs578_de_lfc0_genes, tt642_de_lfc0_genes, file = "R_Data/RNA_PRRX1_DE_Genes.RData")

# Common genes -------------------------------------------------------------------------------------------

hs578_res <- results(hs578_dds, contrast = c("group", "Hs578_TRUE", "Hs578_FALSE"))
tt642_res <- results(tt642_dds, contrast = c("group", "TT642_TRUE", "TT642_FALSE"))
lfc_cutoff <- 1.8

common_genes <- intersect(rownames(hs578_res), rownames(tt642_res))

cbd_res <- data.frame(row.names = common_genes,
                      gene = common_genes,
                      hs578_l2fc = hs578_res[common_genes, "log2FoldChange"],
                      hs578_pval = hs578_res[common_genes, "padj"],
                      tt642_l2fc = tt642_res[common_genes, "log2FoldChange"],
                      tt642_pval = tt642_res[common_genes, "padj"])
cbd_res$label <- with(cbd_res, (abs(hs578_l2fc) > lfc_cutoff) & (abs(tt642_l2fc) > lfc_cutoff) & (tt642_l2fc * hs578_l2fc) > 0)

min_val <- min(c(min(cbd_res$tt642_l2fc), min(cbd_res$hs578_l2fc)))
max_val <- max(c(max(cbd_res$tt642_l2fc), max(cbd_res$hs578_l2fc)))

pdf(file = "analysis/rna_seq_prrx1_kd/summary/diff_genes.pdf")
ggplot(data = cbd_res, aes(x = hs578_l2fc, y = tt642_l2fc, label = gene, color = label)) + 
  geom_point() +
  xlim(c(min_val, max_val)) +
  ylim(c(min_val, max_val)) +
  geom_text_repel(data = subset(cbd_res, label), show.legend = FALSE, max.overlaps = Inf)
dev.off()