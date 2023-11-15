setwd("../")

library(DESeq2)

obj.chip <- load("R_Data/ChIP_Data.RData")
obj.expr <- load("R_Data/RNA_Data.RData")
obj.diff <- load("R_Data/Diff_superenhancers.RData")
obj.heat <- load("R_Data/Heatmap_Metadata.RData")

data <- sweep(chip_se_all_rpkm_norm_uf, 1, rowMeans(chip_se_all_rpkm_norm_uf))
oo <- order(chip_metadata$tnbc.subtype, chip_metadata$cell.line)
data <- data[,oo]
mm <- chip_metadata[oo,]

se_bas_up <- subset(se_basal_diff, log2FoldChange > 0)
se_bas_dn <- subset(se_basal_diff, log2FoldChange < 0)
se_mes_up <- subset(se_mesenchymal_diff, log2FoldChange > 0)
se_mes_dn <- subset(se_mesenchymal_diff, log2FoldChange < 0)
all(rownames(se_bas_up) %in% rownames(data))
all(rownames(se_bas_dn) %in% rownames(data))
all(rownames(se_mes_up) %in% rownames(data))
all(rownames(se_mes_dn) %in% rownames(data))

pdf("analysis/chip_seq/outliers/dse_bar_plots.pdf")
par(mfrow = c(2, 1), mar = c(6, 6, 4, 2))
xx <- data[rownames(se_bas_up), ]
xx <- apply(xx, 2, median)
barplot(xx, col = line_colours$type[mm$tnbc.subtype], las = 2, main = "Basal-specific super-enhancer",
        ylab = "")
title(ylab = "Median H3K27ac expression\ndifference from gene average", line = 3)
xx <- data[rownames(se_mes_up), ]
xx <- apply(xx, 2, median)
barplot(xx, col = line_colours$type[mm$tnbc.subtype], las = 2, main = "Mesenchymal-specific super-enhancers")
title(ylab = "Median H3K27ac expression\ndifference from gene average", line = 3)
dev.off()
