#' Differential expression t-test
ttest_comparison <- function(data, data_ann, group1, group2, path, comparison){
  data1 <- data[,which(colnames(data) %in% group1)]
  data2 <- data[,which(colnames(data) %in% group2)]
  
  data2_mean <- apply(data2, 1, mean)
  data1_mean <- apply(data1, 1, mean)
  data_ann$change <- data2_mean - data1_mean
  
  print("calculating p-values")
  data_ann$p.value <- sapply(1:nrow(data), function(x) t.test(as.numeric(data1[x,]), as.numeric(data2[x,]))$p.value)
  data_ann$padj <- p.adjust(data_ann$p.value, method = "fdr")
  
  sub <- subset(data_ann, padj < 0.05)
  sub <- sub[order(sub$padj),]
  
  write.table(sub, file = paste0(path, "/", comparison, ".csv"),
              sep = ",", row.names = FALSE)
  
  return(sub)
}

# Assess change without stat test to avoid double-dipping
assess_change <- function(data, data_ann, group1, group2, path, comparison){
  data1 <- data[,which(colnames(data) %in% group1)]
  data2 <- data[,which(colnames(data) %in% group2)]
  
  data2_mean <- apply(data2, 1, mean)
  data1_mean <- apply(data1, 1, mean)
  data_ann$change <- data2_mean - data1_mean
  
  sub <- data_ann
  sub <- sub[order(sub$change, decreasing = TRUE), ]
  
  write.table(sub, file = paste0(path, "/", comparison, ".csv"),
              sep = ",", row.names = FALSE)
  
  return(sub)
}