library(DESeq2)

get_deseq_norm_counts <- function(in_counts, in_col_data){
  dds <- DESeqDataSetFromMatrix(countData = in_counts,
                                colData = in_col_data,
                                design = ~ 1)
  keep <- rowSums(counts(dds)) >= 2
  dds <- dds[keep, ]
  dds <- estimateSizeFactors(dds)
  return(counts(dds, normalized = TRUE))
}

setwd("../")

load("R_Data/RNA_PRRX1_OE_Data.RData")
out_path <- "analysis/rna_seq_prrx1_oe/summary/"

cell_lines <- unique(rna_oe_metadata$parental_cell_line)
col_data <- rna_oe_metadata
col_data$group <- paste(col_data$parental_cell_line, col_data$mutant, col_data$condition, sep = "_")

# Counts
for(i in 1:length(cell_lines)){
  cell_line <- cell_lines[i]
  message(cell_line)
  
  # Wild-type counts
  wt_keep <- (col_data$parental_cell_line == cell_line) & (col_data$mutant == "wt")
  wt_counts <- rna_oe_counts[, wt_keep]
  wt_col_data <- col_data[wt_keep, ]
  
  wt_norm_counts <- get_deseq_norm_counts(wt_counts, wt_col_data)
  wt_norm_counts_fname <- paste0(cell_line, "_wt_rna_oe_norm_counts.csv")
  write.table(wt_norm_counts, file = file.path(out_path, wt_norm_counts_fname), 
              sep = ",", col.names = NA)
  
  wt_classes_fname <- paste0(cell_line, "_wt_classes.csv")
  write.table(wt_col_data$group, file.path(out_path, wt_classes_fname), 
              sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Mutant counts
  mut_keep <- (col_data$parental_cell_line == cell_line) & (col_data$mutant == "dH3")
  mut_counts <- rna_oe_counts[, mut_keep]
  mut_col_data <- col_data[mut_keep, ]
  
  mut_norm_counts <- get_deseq_norm_counts(mut_counts, mut_col_data)
  mut_norm_counts_fname <- paste0(cell_line, "_mut_rna_oe_norm_counts.csv")
  write.table(mut_norm_counts, file = file.path(out_path, mut_norm_counts_fname),
              sep = ",", col.names = NA)
  
  mut_classes_fname <- paste0(cell_line, "_mut_classes.csv")
  write.table(mut_col_data$group, file.path(out_path, mut_classes_fname),
              sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Dox counts
  dox_keep <- (col_data$parental_cell_line == cell_line) & (col_data$condition == "dox")
  dox_counts <- rna_oe_counts[, dox_keep]
  dox_col_data <- col_data[dox_keep, ]
  
  dox_norm_counts <- get_deseq_norm_counts(dox_counts, dox_col_data)
  dox_norm_counts_fname <- paste0(cell_line, "_dox_rna_oe_norm_counts.csv")
  write.table(dox_norm_counts, file = file.path(out_path, dox_norm_counts_fname),
              sep = ",", col.names = NA)
  
  dox_classes_fname <- paste0(cell_line, "_dox_classes.csv")
  write.table(dox_col_data$group, file.path(out_path, dox_classes_fname),
              sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Metadata
write.table(rna_oe_metadata, file = file.path(out_path, "rna_oe_metadata.csv"), sep = ",", row.names = FALSE)

# RPKM norm
write.table(rna_oe_rpkm_norm_uf, file = file.path(out_path, "rna_oe_fpkm_norm_unfiltered.csv"), sep = ",", col.names = NA)
