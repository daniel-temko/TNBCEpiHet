source("../utilities/utils.R")
setwd("../")
library(reshape2)

extract_vals <- function(df, row_ids, 
                         col_ids, 
                         row_names = c("high", "low"), 
                         col_names = c("bas", "lum", "mes")){
  df <- df[row_ids, col_ids]
  rownames(df) <- row_names
  colnames(df) <- col_names
  return(df)
}

sum_df <- read.table("analysis/integration/summary/Summary_Diff_DE_genes_RNA_ChIP.csv", 
                        fill = TRUE,
                        header = TRUE,
                        sep = ",")

mat_list <- list(de_genes = extract_vals(sum_df, 1:2, 2:4),
                 de_tfs = extract_vals(sum_df, 5:6, 2:4),
                 dse_genes = extract_vals(sum_df, 1:2, 7:9),
                 dse_tfs = extract_vals(sum_df, 5:6, 7:9))

ids <- c("de_genes_tab",
         "de_tfs_tab",
         "dse_genes_tab",
         "dse_tfs_tab")

p_vals <- sapply(mat_list, function(x) apply(x, 2, function(y) {
  binom.test(y, p=0.5)$p.value
}))

p_vals_lg <- melt(p_vals)
p_vals_lg$p.adj.holm <- p.adjust(p_vals_lg$value, method="holm")
p_vals_lg <- p_vals_lg[order(p_vals_lg$p.adj.holm),]
colnames(p_vals_lg)[1:3] <- c("type", "comparison", "p.value")

write_csv(p_vals_lg, file = "analysis/integration/summary/binom_test_p_vals.csv")
