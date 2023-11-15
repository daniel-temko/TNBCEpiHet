source("../utilities/integration.R")
library(ggplot2)

setwd("../")

add_transparency <- function(hex_col, alpha = 255){
  do.call(rgb, as.list(c(col2rgb(hex_col), alpha, maxColorValue = 255)))
}

load("R_Data/RNA_Data.RData")
load("R_Data/Methyl_Data.RData")
load("R_Data/ChIP_Data.RData")
load("R_Data/Gene_SE_Data.RData")

comp_levels <- c("rna_vs_methyl_se",
                 "rna_vs_methyl_gb",
                 "rna_vs_methyl_tss",
                 "rna_vs_chip_se")
cols <- c(rna_vs_methyl_se = "#7C00FF", 
          rna_vs_methyl_gb = add_transparency("#7C00FF00", 205),
          rna_vs_methyl_tss = add_transparency("#7C00FF00", 127),
          rna_vs_chip_se = "#FF8000")

# Main analysis with unfiltered RNA data

# RNA and DNA methylation in intergenic SE's
rms_df <- align_and_correlate(rna_rpkm_norm_uf, methyl_gene_se_agg_m, "rna_vs_methyl_se")

# RNA and DNA methylation in gene bodies
rmg_df <- align_and_correlate(rna_rpkm_norm_uf, methyl_gb_agg_m, "rna_vs_methyl_gb")

# RNA and DNA methylation in promoters
rmt_df <- align_and_correlate(rna_rpkm_norm_uf, methyl_tss_agg_m, "rna_vs_methyl_tss")

# RNA and ChIP-seq
rcgs_df <- align_and_correlate(rna_rpkm_norm_uf, chip_gene_se_rpkm_norm_uf, "rna_vs_chip_se")

res_df <- rbind(rms_df, rmg_df, rmt_df, rcgs_df)

med_df <- aggregate(cor ~ comp, res_df, median)

p_val_df <- data.frame(p_val = sapply(comp_levels, function(x) {
  sub <- subset(res_df, comp == x)
  wilcox.test(sub$cor)$p.value
}))
p_val_df$median <- sapply(comp_levels, function(x) {
  sub <- subset(res_df, comp == x)
  median(sub$cor)
})
p_val_df <- p_val_df[,c("median", "p_val")]

write.table(p_val_df, 
            "analysis/integration/overall_regulation/rna_methyl_chip_p_vals.csv",
            col.names = NA,
            sep = ",")

plot_df <- res_df
plot_df$comp <- factor(plot_df$comp, levels = comp_levels)
plot_ann <- data.frame(x = 1:length(comp_levels),
                       y = 1,
                       comp = comp_levels, 
                       n = paste0("N=", format(as.numeric(table(plot_df$comp)), big.mark = ",")))

pdf("analysis/integration/overall_regulation/rna_methyl_chip.pdf")
ggplot(plot_df, aes(x = comp, y = cor, fill = comp)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1,
               show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Comparison") + 
  ylab("Pearson Correlation") +
  geom_text(size = 6,
            data = plot_ann,
            aes(x = x, y = y, label = n)) +
  scale_fill_manual(values = cols)
dev.off()

# Sensitivity using filtered RNA data

# RNA and DNA methylation in intergenic SE's
rms_df2 <- align_and_correlate(rna_rpkm_norm, methyl_gene_se_agg_m, "rna_vs_methyl_se")

# RNA and DNA methylation in gene bodies
rmg_df2 <- align_and_correlate(rna_rpkm_norm, methyl_gb_agg_m, "rna_vs_methyl_gb")

# RNA and DNA methylation in promoters
rmt_df2 <- align_and_correlate(rna_rpkm_norm, methyl_tss_agg_m, "rna_vs_methyl_tss")

# RNA and ChIP-seq
rcgs_df2 <- align_and_correlate(rna_rpkm_norm, chip_gene_se_rpkm_norm_uf, "rna_vs_chip_se")

res_df2 <- rbind(rms_df2, rmg_df2, rmt_df2, rcgs_df2)

p_val_df2 <- data.frame(p_val = sapply(comp_levels, function(x){
  sub <- subset(res_df2, comp == x)
  wilcox.test(sub$cor)$p.value
}))
p_val_df2$median <- sapply(comp_levels, function(x){
  sub <- subset(res_df2, comp == x)
  median(sub$cor)
})
p_val_df2 <- p_val_df2[,c("median", "p_val")]

write.table(p_val_df2,
            "analysis/integration/overall_regulation/rna_methyl_chip_p_vals_filtered.csv",
            col.names = NA,
            sep = ",")

plot_df2 <- res_df2
plot_df2$comp <- factor(plot_df2$comp, levels = comp_levels)
plot_ann2 <- data.frame(x = 1:length(comp_levels),
                        y = 1,
                        comp = comp_levels,
                        n = paste0("N=", format(as.numeric(table(plot_df2$comp)), big.mark = ",")))

pdf("analysis/integration/overall_regulation/rna_methyl_chip_filtered.pdf")
p <- ggplot(plot_df2, aes(x = comp, y = cor, fill = comp)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1,
               show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Comparison") +
  ylab("Pearson Correlation") +
  geom_text(size = 6,
            data = plot_ann2,
            aes(x = x, y = y, label = n)) +
  scale_fill_manual(values = cols)
print(p)
dev.off()
