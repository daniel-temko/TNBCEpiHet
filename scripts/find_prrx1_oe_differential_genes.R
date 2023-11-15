library(DESeq2)
library(ggplot2)
library(ggrepel)

setwd("../")

# Load data -------------------------------------------------------------------------------------------

load("R_Data/RNA_PRRX1_OE_Data.RData")
out_path <- "analysis/rna_seq_prrx1_oe/differential_genes/"
cell_lines <- unique(rna_oe_metadata$parental_cell_line)

contrasts <- data.frame(a = c(paste0(cell_lines, "_wt_dox"),
                              paste0(cell_lines, "_dH3_dox"),
                              paste0(cell_lines, "_wt_dox")),
                        b = c(paste0(cell_lines, "_wt_no_dox"),
                              paste0(cell_lines, "_dH3_no_dox"),
                              paste0(cell_lines, "_dH3_dox")))

# DESeq2 -------------------------------------------------------------------------------------------

col_data <- rna_oe_metadata
col_data$group <- paste(col_data$parental_cell_line, col_data$mutant, col_data$condition, sep = "_")
col_data$group <- factor(col_data$group)

for(i in 1:nrow(contrasts)){
  comp_id <- paste0(contrasts[i,], collapse = "_vs_")
  print(comp_id)
  
  col_data_sub <- droplevels(subset(col_data, group %in% contrasts[i,]))  
  counts_sub <- rna_oe_counts[, rownames(col_data_sub)]
  
  dds <- DESeqDataSetFromMatrix(countData = counts_sub,
                                colData = col_data_sub,
                                design =  ~ group)
  keep <- rowSums(counts(dds)) > 0
  dds <- dds[keep, ]
  
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("group", contrasts$a[i], contrasts$b[i]))
  
  sub <- subset(res, padj < 0.05)
  sub <- sub[order(sub$padj), ]
  
  write.table(sub, 
              file = paste0(out_path, comp_id, "_lfc0_genes.csv"),
              sep = ",",
              col.names = NA)
  
  assign(paste0("diff_", comp_id), res)
}

comp_ids <- apply(contrasts, 1, function(x) paste0(x, collapse = "_vs_"))
save(list = paste0("diff_", comp_ids), file = "R_Data/RNA_PRRX1_OE_DE_Genes.RData")

