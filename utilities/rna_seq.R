library(DESeq2)

#' Takes count_data, and sample_table and returns differentially expressed genes in 
#' foreground vs. background condition in terms of the condition column of sample_table
#' 
#' Note: Foreground and background conditions should be strings
#' Note: condition column of sample_table should be a factor - order of levels does not matter
#' Note: Genes with positive log2fc are higher in foreground condition
de_genes <- function(count_data, count_data_ann, sample_table, foreground_condition, background_condition, fdr=NULL){
  # Subset count_data and sample_table to conditions to be tested
  keep <- sample_table$condition %in% c(foreground_condition, background_condition)
  sample_table <- droplevels(sample_table[keep, ])
  count_data <- count_data[, keep]
  
  message("running with ", ncol(count_data), " samples")
  
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = sample_table,
                                design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", foreground_condition, background_condition))
  res <- as.data.frame(res)
  
  # Add gene annotations
  res$gene.symbol <- rownames(res)
  res <- cbind(res, count_data_ann[rownames(res), , drop = FALSE])
  
  # Reorder columns
  res <- bring_multiple_to_front(res, c("gene.symbol", colnames(count_data_ann)))
  
  # Sort row by padj and filter
  res <- res[order(res$padj), ]
  if(!is.null(fdr)){
    message("filtering to padj<", fdr)
    res <- subset(res, padj < 0.05)  
  } else {
    message("no fdr value specified, not filtering")
  }
  
  return(res)
}
