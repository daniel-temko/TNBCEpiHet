setwd("../")

source("scripts/utilities/supervised_analysis.R")

load("R_Data/Metab_Data.RData")

tnbc_types <- unique(metab_rep_metadata$tnbc.subtype)

for(i in 1:length(tnbc_types)){
  type <- tnbc_types[i]
  print(type)
  group1 <- rownames(metab_rep_metadata)[which(metab_rep_metadata$tnbc.subtype != type)]
  group2 <- rownames(metab_rep_metadata)[which(metab_rep_metadata$tnbc.subtype == type)]
  
  file_name <- paste0(type, "_de_metabolites")
  
  stopifnot(all(rownames(metab_ann) == rownames(metab_rep_norm)))
  ttest_comparison(metab_rep_norm, 
                   metab_ann, 
                   group1, 
                   group2, 
                   "analysis/metabolomics/differential_metabolites", 
                   file_name)
}

