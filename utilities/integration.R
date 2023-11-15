#' Filter and reorder the rows and columns of data1 and data2 to match each 
#' other. Remove any rows that have zero variance in either data-set. 
#' Return the row-wise correlations for the filtered, ordered data-sets
align_and_correlate <- function(data1, data2, comp_id){
  #' Filter and reorder. Coerce to matrices to enable taking row-wise
  #' correlations
  common_genes <- sort(intersect(rownames(data1), rownames(data2)))
  common_cols <- sort(intersect(colnames(data1), colnames(data2)))
  data1 <- as.matrix(data1[common_genes, common_cols])
  data2 <- as.matrix(data2[common_genes, common_cols])
  keep1 <- apply(data1, 1, var) > 0
  keep2 <- apply(data2, 1, var) > 0
  keep <- keep1 & keep2
  data1 <- data1[keep, ]
  data2 <- data2[keep, ]
  common_genes <- common_genes[keep]
  # Calculate correlations per row
  cor_vals <- sapply(1:length(common_genes), function(x){
    cor(data1[x, ], data2[x, ])
  })
  data.frame(gene = common_genes, cor = cor_vals, comp = comp_id)
}