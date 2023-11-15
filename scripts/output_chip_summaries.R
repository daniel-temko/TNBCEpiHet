library(reshape2)

setwd("../")

load("R_Data/ChIP_Data.RData")
load("R_Data/ChIP_Overlap_Data.RData")
load("R_Data/Heatmap_Metadata.RData")

canon <- paste0("chr", c(1:22, "X", "Y"))
ordered_types <- c("mesenchymal", "luminal", "basal")

all(rownames(chip_se_olp_peaks) == rownames(chip_se_all_rpkm_norm_uf))
all(rownames(chip_se_ann) == rownames(chip_se_all_rpkm_norm_uf))
all(rownames(chip_se_olp_peaks) == rownames(chip_se_all_counts))

# Get ranked lists of SE's present (overlapped in) each cell line
ranked_se_list <- lapply(as.character(chip_metadata$cell.line), function(x){
  df <- data.frame(row.names = rownames(chip_se_all_counts),
                   cell.line = x,
                   counts = chip_se_all_counts[[x]],
                   chr.region = chip_se_ann$chr.region,
                   gene.symbol = chip_se_ann$gene.symbol,
                   all.gene.symbols = chip_se_olp_peaks2$all.closest.genes)
  keep <- (chip_se_olp_peaks[[x]] > 0)
  df <- df[keep, ]
  keep <- df$chr.region %in% canon
  df <- df[keep, ]
  df <- df[order(df$counts, decreasing = TRUE), ]
  df$se.rank <- 1:nrow(df)
  df
})
names(ranked_se_list) <- as.character(chip_metadata$cell.line)

ranked_se_data <- do.call(rbind, ranked_se_list)
ranked_se_data$se.rank <- factor(ranked_se_data$se.rank, levels = 1:max(ranked_se_data$se.rank))
ranked_se_data$cell.line <- factor(ranked_se_data$cell.line, levels = chip_metadata$cell.line)

se_ranks <- dcast(ranked_se_data, se.rank ~ cell.line, value.var = "gene.symbol", fill = "")

se_ranks2 <- dcast(ranked_se_data, se.rank ~ cell.line, value.var = "all.gene.symbols", fill = "")

save(se_ranks, se_ranks2, file = "R_Data/ChIP_SE_Ranks.RData")

write.table(se_ranks, file = "analysis/chip_seq/summary/tnbc_se_summary.csv", sep = ",", row.names = FALSE)

write.table(se_ranks2, file = "analysis/chip_seq/summary/tnbc_se_summary_all_genes.csv", sep = ",", row.names = FALSE)
