source("../utilities/plots.R")
source("../utilities/utils.R")
library(FSA)

setwd("../")

# Load data -----------------------------------------------------------------------------------------

load("R_Data/ChIP_Data.RData")
load("R_Data/Heatmap_Metadata.RData")

chip_titles <- data.frame(title = c("Peaks", "Superenhancers", "Enhancers"),
                          type = c("peak", "superenhancer", "enhancer"),
                          type2 = c("peaks", "superenhancers", "enhancers"),
                          id = c("peak", "se", "enh"))

# Analysis -----------------------------------------------------------------------------------------

n_read_vals <- list(peak = colSums(chip_peak_all_counts),
                    se = colSums(chip_se_all_counts),
                    enh = colSums(chip_enh_all_counts))
prop_read_vals <- lapply(n_read_vals, function(x) {
  x / chip_mapped_reads
})

prop_read_aov_list <- lapply(prop_read_vals, function(x){
  aov(x ~ chip_metadata$tnbc.subtype)
})

(prop_read_p <- sapply(prop_read_aov_list, function(x) summary(x)[[1]][1, 5]))

prop_read_kw_list <- lapply(prop_read_vals, function(x){
  kruskal.test(x ~ chip_metadata$tnbc.subtype)
})

(prop_read_kw_p <- sapply(prop_read_kw_list, function(x) x$p.value))

rpkm_norm_vals <- list(peak = colMeans(chip_peak_all_rpkm_norm_uf),
                       se = colMeans(chip_se_all_rpkm_norm_uf),
                       enh = colMeans(chip_enh_all_rpkm_norm_uf))

rpkm_norm_aov_list <- lapply(rpkm_norm_vals, function(x) {
  aov(x ~ chip_metadata$tnbc.subtype)
})

(rpkm_norm_p <- sapply(rpkm_norm_aov_list, function(x) summary(x)[[1]][1, 5]))

rpkm_norm_kw_list <- lapply(rpkm_norm_vals, function(x) {
  kruskal.test(x ~ chip_metadata$tnbc.subtype)
})

(rpkm_norm_kw_p <- sapply(rpkm_norm_kw_list, function(x) x$p.value))


# Plots -----------------------------------------------------------------------------------------

for(i in 1:length(rpkm_norm_vals)){
  ylab <- paste0("Proportion of H3K27ac Reads in ", chip_titles$title[i])
  
  pdf(paste0("analysis/chip_seq/total_signal/", chip_titles$type[i], "_h3k27ac_prop_reads_by_group.pdf"), width = 9, height = 9)
  plot_by_group(prop_read_vals[[i]], chip_metadata$tnbc.subtype, ylab)
  dev.off()
  
  pdf(paste0("analysis/chip_seq/total_signal/", chip_titles$type[i], "_h3k27ac_prop_reads_by_line.pdf"))
  plot_by_line(prop_read_vals[[i]], chip_metadata$tnbc.subtype, ylab, chip_metadata$cell.line)
  dev.off()
}

for(i in 1:length(rpkm_norm_vals)){
  ylab <- paste0("Average H3K27ac Signal in ", chip_titles$title[i])
  
  pdf(paste0("analysis/chip_seq/total_signal/", chip_titles$type[i], "_h3k27ac_by_group.pdf"), width = 9, height = 9)
    plot_by_group(rpkm_norm_vals[[i]], chip_metadata$tnbc.subtype, ylab)
  dev.off()
  
  pdf(paste0("analysis/chip_seq/total_signal/", chip_titles$type[i], "_h3k27ac_by_line.pdf"))
    plot_by_line(rpkm_norm_vals[[i]], chip_metadata$tnbc.subtype, ylab, chip_metadata$cell.line)
  dev.off()
}

# Save results -----------------------------------------------------------------------------------------

prop_read_p_vals <- data.frame(test = chip_titles$id,
                                p.val = prop_read_kw_p)

write_csv(prop_read_p_vals, "analysis/chip_seq/total_signal/chip_total_prop_read_p_vals.csv")

rpkm_norm_p_vals <- data.frame(test = chip_titles$id,
                               p.val = rpkm_norm_kw_p)

write_csv(rpkm_norm_p_vals, "analysis/chip_seq/total_signal/chip_total_p_vals.csv")

for(i in 1:length(prop_read_kw_list)){
  prop_read_dunn_test <- dunnTest(prop_read_vals[[i]] ~ chip_metadata$tnbc.subtype, method = "holm")
  fn <- paste0(chip_titles$type2[i], "_prop_reads_pairwise_comparisons.csv")
  write.table(prop_read_dunn_test[[2]], 
              paste0("analysis/chip_seq/total_signal/", fn),
              sep = ",",
              col.names = NA)
}

for(i in 1:length(rpkm_norm_kw_list)){
  rpkm_norm_dunn_test <- dunnTest(rpkm_norm_vals[[i]] ~ chip_metadata$tnbc.subtype, method = "holm")
  fn <- paste0(chip_titles$type2[i], "_pairwise_comparisons.csv")
  write.table(rpkm_norm_dunn_test[[2]], 
              paste0("analysis/chip_seq/total_signal/", fn),
              sep = ",",
              col.names = NA)
}
