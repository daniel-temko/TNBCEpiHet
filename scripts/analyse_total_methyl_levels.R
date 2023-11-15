source("../utilities/utils.R")
source("../utilities/plots.R")

setwd("../")

# Load data -----------------------------------------------------------------------------------------

load("R_Data/Methyl_Data.RData")
load("R_Data/Heatmap_Metadata.RData")

m_titles <- data.frame(title = c("All", "Int SE", "Superenhancers", "Enhancers", "Gene Bodies", "TSS"),
                       type = c("all", "superenhancer_int", "superenhancer", "enhancer", "gene_body", "tss"),
                       id = c("all", "se_int", "se", "enh", "gb", "tss"))

# Analysis -----------------------------------------------------------------------------------------

m_vals <- list(all = colMeans(methyl_all_m),
                     se.int = colMeans(methyl_se_int_agg_m),
                     se = colMeans(methyl_se_agg_m),
                     enh = colMeans(methyl_enh_agg_m),
                     gb = colMeans(methyl_gb_agg_m),
                     tss = colMeans(methyl_tss_agg_m))

m_aov_list <- lapply(m_vals, function(x) {
  aov(x ~ methyl_metadata$tnbc.subtype)
})

m_p <- sapply(m_aov_list, function(x) summary(x)[[1]][1,5])

m_kw_list <- lapply(m_vals, function(x) {
  kruskal.test(x ~ methyl_metadata$tnbc.subtype)
})

m_kw_p <- sapply(m_kw_list, function(x) x$p.value)

# Plots -----------------------------------------------------------------------------------------

for(i in 1:length(m_vals)){
  ylab <- paste0("Average M-value in ", m_titles$title[i])
  cell_lines <- methyl_metadata$Cell.Line
  
  pdf(paste0("analysis/dna_methylation/total_signal/", m_titles$type[i], "_m_vals_by_group.pdf"), width = 9, height = 9)
    plot_by_group(m_vals[[i]], methyl_metadata$tnbc.subtype, ylab)
  dev.off()
  
  pdf(paste0("analysis/dna_methylation/total_signal/", m_titles$type[i], "_m_vals_by_line.pdf"))
    plot_by_line(m_vals[[i]], methyl_metadata$tnbc.subtype, ylab, cell_lines)
  dev.off()
}

# Save results -----------------------------------------------------------------------------------------

p_vals <- data.frame(test = m_titles$id,
                     p.val = m_kw_p)
write_csv(p_vals, "analysis/dna_methylation/total_signal/m_val_p_vals.csv")
