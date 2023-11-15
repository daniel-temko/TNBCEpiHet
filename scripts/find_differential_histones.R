setwd("../")

source("scripts/utilities/supervised_analysis.R")
source("scripts/utilities/preprocessing.R")

obj.hist <- load("R_Data/HMS_Data.RData")
obj.clst <- load("analysis/hms/exploratory/r_data/raw_data_euc_top20pc_clustering.RData") # col_hclust

tnbc_types <- unique(hms_metadata$tnbc.subtype)

for(i in 1:length(tnbc_types)){
  type <- tnbc_types[i]
  print(type)
  group1 <- rownames(hms_metadata)[which(hms_metadata$tnbc.subtype != type)]
  group2 <- rownames(hms_metadata)[which(hms_metadata$tnbc.subtype == type)]
  
  file_name <- paste0(type, "_diff_marks")
  
  stopifnot(all(rownames(hms_norm) == rownames(hms_ann)))
  ttest_comparison(hms_norm, 
                   hms_ann, 
                   group1, 
                   group2, 
                   "analysis/hms/differential_marks", 
                   file_name)
}

#' Compare two major clusters
cluster_labels_k2 <- cutree(col_hclust, k = 2)
stopifnot(all(names(cluster_labels_k2) == colnames(hms_norm)))
# Switch labels so that 2 refers to right-hand cluster
#' left: "1", right: "2"
cc <- 3 - cluster_labels_k2

hms_variances <- variable_features(hms_norm)
n_feat <- round(nrow(hms_norm) * 0.2)
xx <- hms_norm[hms_variances[1:n_feat],]
xx <- sweep(xx, 1, rowMeans(xx))
aa <- hms_ann[hms_variances[1:n_feat], , drop = FALSE]

group1 <- names(cc)[which(cc == 1)]
group2 <- names(cc)[which(cc == 2)]

assess_change(xx,
              aa,
              group1,
              group2,
              "analysis/hms/cluster_comparison",
              "hms_cluster_2vs1_mark_changes")

for(i in 1:nrow(xx)){
  mm <- rownames(xx)[i]
  nn <- gsub("\\)", "_", gsub("\\(", "_", mm))
  pdf(paste0("analysis/hms/cluster_comparison/", i, "_", nn, ".pdf"))
  boxplot(as.numeric(xx[i,]) ~ cc, main = mm, ylab = mm, xlab = "Cluster")
  dev.off()
}

