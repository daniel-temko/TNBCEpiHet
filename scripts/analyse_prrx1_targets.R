source("../utilities/chip_seq.R")
library(DESeq2)
library(UpSetR)

setwd("../")

# Load data ---------------------------------------------------------------------------------------------------------------------------------------

load(file = "R_Data/PRRX1_ChIP_Data.RData")
load(file = "R_Data/RNA_PRRX1_DE_Genes.RData")

hs578_all <- as.data.frame(hs578_de_lfc0_genes)
tt642_all <- as.data.frame(tt642_de_lfc0_genes)

# Create target summaries ---------------------------------------------------------------------------------------------------

prrx1_chip_df <- annotate_targets(prrx1_chip_targets, prrx1_chip_metadata$cell.line)
prrx1_chip_df <- add_expression(prrx1_chip_df, hs578_all, "hs578_kd")
prrx1_chip_df <- add_expression(prrx1_chip_df, tt642_all, "TT642_kd")
prrx1_chip_cl_df <- annotate_targets(prrx1_chip_cl_targets, colnames(prrx1_chip_cl_targets))
prrx1_chip_cl_df <- add_expression(prrx1_chip_cl_df, hs578_all, "hs578_kd")
prrx1_chip_cl_df <- add_expression(prrx1_chip_cl_df, tt642_all, "TT642_kd")

df <- prrx1_chip_df
colnames(df)[match("num_lines", colnames(df))] <- "num_targets"

hs578_comb_all <- cbind(hs578_all, df[rownames(hs578_all), c("targets", "num_targets")])
hs578_na_rows <- which(is.na(hs578_comb_all$targets))
hs578_comb_all[hs578_na_rows, c("targets", "num_targets")] <- 0
hs578_comb_all <- hs578_comb_all[order(hs578_comb_all$num_targets, decreasing = TRUE),]

tt642_comb_all <-cbind(tt642_all, df[rownames(tt642_all), c("targets", "num_targets")])
tt642_na_rows <- which(is.na(tt642_comb_all$targets))
tt642_comb_all[tt642_na_rows, c("targets", "num_targets")] <- 0
tt642_comb_all <- tt642_comb_all[order(tt642_comb_all$num_targets, decreasing = TRUE),]

write.table(hs578_comb_all, "analysis/prrx1_chip_seq/targets/hs578_lfc0_with_targets.csv", sep = ",", col.names = NA)
write.table(tt642_comb_all, "analysis/prrx1_chip_seq/targets/tt642_lfc0_with_targets.csv", sep = ",", col.names = NA)
write.table(prrx1_chip_df, "analysis/prrx1_chip_seq/targets/prrx1_chip_target.csv", sep = ",", col.names = NA)
write.table(prrx1_chip_cl_df, "analysis/prrx1_chip_seq/targets/prrx1_chip_cl_target.csv", sep = ",", col.names = NA)

# Upset plots ---------------------------------------------------------------------------------------------------

temp <- list(hs578.target = rownames(prrx1_chip_cl_df)[prrx1_chip_cl_df$Hs578],
             hs578.kd.up = rownames(hs578_all)[hs578_all$log2FoldChange > 0],
             hs578.kd.dn = rownames(hs578_all)[hs578_all$log2FoldChange < 0])
pdf("analysis/prrx1_chip_seq/targets/hs578_chip_kd_overlap.pdf")
  upset(fromList(temp), 
        order.by = "freq", 
        point.size = 3.5, 
        line.size = 2, 
        sets.x.label = "Number of Genes",
        mainbar.y.label = "Gene Interactions",
        text.scale = c(1.3, 1.3, 1, 1, 2, 1.5))
dev.off()

temp <- list(TT642.target = rownames(prrx1_chip_cl_df)[prrx1_chip_cl_df$TT642],
             TT642.kd.up = rownames(tt642_all)[tt642_all$log2FoldChange > 0],
             TT642.kd.dn = rownames(tt642_all)[tt642_all$log2FoldChange < 0])
pdf("analysis/prrx1_chip_seq/targets/tt642_chip_kd_overlap.pdf")
upset(fromList(temp), 
      order.by = "freq", 
      point.size = 3.5, 
      line.size = 2, 
      sets.x.label = "Number of Genes",
      mainbar.y.label = "Gene Interactions",
      text.scale = c(1.3, 1.3, 1, 1, 2, 1.5))
dev.off()
