source("../utilities/utils.R")
source("../utilities/rna_seq.R")

setwd("../")

load("R_Data/RNA_Data.RData") # rna_counts, rna_counts_ann
load("analysis/hms/exploratory/r_data/raw_data_euc_top20pc_clustering.RData") # col_hclust

#' Compare two major clusters
#' left: "2", right: "1"
cluster_labels_k2 <- cutree(col_hclust, k = 2)

sample_table <- data.frame(sample = colnames(rna_counts),
                           condition = factor(cluster_labels_k2[colnames(rna_counts)]))

#' Note: The right-hand cluster number is "1", thus "1" is used as foreground condition, 
#' however everywhere else (i.e. in 2vs1) 2 refers to the right-hand cluster by convention
res_2vs1 <- de_genes(rna_counts, rna_counts_ann[, "dbd.tf", drop = FALSE], sample_table, "1", "2", 0.05)

fname_2vs1 <- "analysis/rna_seq/differential_genes/hms_cluster_2vs1_de_lfc0_genes.csv"
write.table(res_2vs1, file = fname_2vs1, sep = ",", row.names = FALSE)

#' Compare two sub-clusters to each-other within each of the two major clusters
#' far-left: 4, near-left: 3, near-right: 1, far-right: 2
cluster_labels_k4 <- cutree(col_hclust, k = 4)

sample_table <- data.frame(sample = colnames(rna_counts),
                           condition = factor(cluster_labels_k4[colnames(rna_counts)]))

res_1.2vs1.1 <- de_genes(rna_counts, rna_counts_ann[, "dbd.tf", drop = FALSE], sample_table, "3", "4", 0.05)

fname_1.2vs1.1 <- "analysis/rna_seq/differential_genes/hms_cluster_1.2vs1.1_de_lfc0_genes.csv"
write.table(res_1.2vs1.1, file = fname_1.2vs1.1, sep = ",", row.names = FALSE)

res_2.2vs2.1 <- de_genes(rna_counts, rna_counts_ann[, "dbd.tf", drop = FALSE], sample_table, "2", "1", 0.05)

fname_2.2vs2.1 <- "analysis/rna_seq/differential_genes/hms_cluster_2.2vs2.1_de_lfc0_genes.csv"
write.table(res_2.2vs2.1, file = fname_2.2vs2.1, sep = ",", row.names = FALSE)