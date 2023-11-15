setwd("../")

source("scripts/utilities/supervised_analysis.R")

call_dmr <- function(m_values, ann, long_name, sample_table, design, path = "analysis/dna_methylation/dmr/"){
  for(i in 1:nrow(design)){
    print(i)
    
    group1 <- colnames(m_values)[which(sample_table$tnbc.subtype != design$group[i])]
    group2 <- colnames(m_values)[which(sample_table$tnbc.subtype == design$group[i])]
    
    comparison_name <- paste0(long_name, "_", design$group[i], "_dmrs")
    
    stopifnot(all(rownames(ann) == rownames(m_values)))
    ttest_comparison(m_values, ann, group1, group2, path, comparison_name)
  }
}

# Find DMRs ----------------------------------------------------------------------------------

load("R_Data/Methyl_Data.RData")

design <- data.frame(group = c("mesenchymal", "basal", "luminal"))

call_dmr(methyl_se_agg_m, methyl_se_agg_ann, "superenhancers", methyl_metadata, design)
call_dmr(methyl_se_int_agg_m, methyl_se_int_agg_ann, "int_superenhancers", methyl_metadata, design)
call_dmr(methyl_enh_agg_m, methyl_enh_agg_ann, "enhancers", methyl_metadata, design)
call_dmr(methyl_gb_agg_m, methyl_gb_agg_ann, "gene_body", methyl_metadata, design)
call_dmr(methyl_tss_agg_m, methyl_tss_agg_ann, "tss", methyl_metadata, design)


