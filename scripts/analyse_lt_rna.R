library(DESeq2)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(UpSetR)
library(eulerr)

setwd("../")

# Load data -------------------------------------------------------------------------------------------

load("R_Data/RNA_PRRX1_LT_Data.RData")

treatment_lengths <- c("4 weeks", "8 weeks")
cell_lines <- c("Hs578", "TT642")
nt_map <- c(no_sh = "NT-control", shPRRX1_1 = "shPRRX1_1", shPRRX1_3 = "shPRRX1_3")

rna_lt_metadata$shrna2 <- nt_map[rna_lt_metadata$shrna]

de_dir <- "analysis/rna_seq_prrx1_lt/differential_genes/"

# PRRX1 knockdown -------------------------------------------------------------------------------------------

for(i in 1:length(treatment_lengths)){
  lgth <- treatment_lengths[i]
  ids <- which(rna_lt_metadata$treatment.length == lgth)
  expr <- rna_lt_rpkm[,ids]
  expr_norm <- rna_lt_rpkm_norm[,ids]
  samp_tab <- rna_lt_metadata[ids,]
  
  data <- data.frame(prrx1.rpkm = as.numeric(expr["PRRX1",]), 
                     prrx1.log.rpkm = as.numeric(expr_norm["PRRX1",]))
  data$condition <- factor(samp_tab$condition, levels = c("no_DOX", "DOX"))
  data$experiment <- with(samp_tab, paste(parental.cell.line, shrna2))
  
  pdf(paste0("analysis/rna_seq_prrx1_lt/comparison/prrx1_rpkm_", lgth, ".pdf"))
  p <- ggplot(data, aes(x = experiment, y = prrx1.rpkm, fill = condition)) +
    geom_bar(stat = "identity", position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
  dev.off()
  
  pdf(paste0("analysis/rna_seq_prrx1_lt/comparison/prrx1_log_rpkm_", lgth, ".pdf"))
  p <- ggplot(data, aes(x = experiment, y = prrx1.log.rpkm, fill = condition)) +
    geom_bar(stat = "identity", position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
  dev.off()
}

# DESeq2 -------------------------------------------------------------------------------------------

col_data <- rna_lt_metadata
rownames(col_data) <- col_data$sample
col_data$treat <- ! col_data$shrna == "no_sh"
col_data$induced <- col_data$treat & col_data$condition == "DOX"
col_data$group <- factor(col_data$induced)

design <- expand.grid(cell.line = cell_lines, trt.length = treatment_lengths)

res_list <- list()
for(i in 1:nrow(design)){
  cell_line <- design$cell.line[i]
  trt_length <- design$trt.length[i]
  id <- tolower(gsub(" ", "_", paste(cell_line, trt_length)))
  print(id)
  
  col_data_sub <- droplevels(subset(col_data, 
                                    (parental.cell.line == cell_line) & 
                                    (treatment.length == trt_length)))
  counts_sub <- rna_lt_counts[,rownames(col_data_sub)]
  
  dds <- DESeqDataSetFromMatrix(countData = counts_sub,
                                colData = col_data_sub,
                                design = ~ group)
  
  keep <- rowSums(counts(dds)) > 0
  dds <- dds[keep, ]
  
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("group", "TRUE", "FALSE"))
  
  sub <- subset(res, padj < 0.05)
  sub <- sub[order(sub$padj), ]
  
  write.table(sub, 
              file = paste0(de_dir, cell_line, "_", gsub(" ", "_", trt_length), "_de_lfc0_genes.csv"),
              sep = ",",
              col.names = NA)
  res_list[[id]] <- sub
}
save(res_list, file = "R_Data/RNA_PRRX1_LT_DE_Genes.RData")

# Venn diagrams of differential expression overlaps at different time points
common_4 <- fromList(list(hs578_4wk = rownames(res_list$hs578_4_weeks), 
                          tt642_4wk = rownames(res_list$tt642_4_weeks)))
x <- euler(common_4)
pdf("analysis/rna_seq_prrx1_lt/comparison/diff_genes_overlap_4_weeks.pdf")
plot(x, labels = list(fontsize = 20))
dev.off()

common_8 <- fromList(list(hs578_8 = rownames(res_list$hs578_8_weeks),
                          tt642_8wk = rownames(res_list$tt642_8_weeks)))
x <- euler(common_8)
pdf("analysis/rna_seq_prrx1_lt/comparison/diff_genes_overlap_8_weeks.pdf")
plot(x, labels = list(fontsize = 20))
dev.off()

intersect(rownames(res_list$hs578_4_weeks), rownames(res_list$tt642_4_weeks))
intersect(rownames(res_list$hs578_8_weeks), rownames(res_list$tt642_8_weeks))
