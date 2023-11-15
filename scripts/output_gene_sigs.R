setwd("../")

load("R_Data/Gene_Sigs.RData")

# TNBC type sigs
for(i in 1:length(tnbc_gsl)){
  id <- names(tnbc_gsl)[i]
  
  out_up <- data.frame(gene = tnbc_gsl[[i]]$up)
  fname_up <- paste0(id, "_sig_up.csv")
  write.table(out_up, paste0("analysis/integration/summary/", fname_up), 
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              sep = ",")
  
  out_dn <- data.frame(gene = tnbc_gsl[[i]]$dn)
  fname_dn <- paste0(id, "_sig_dn.csv")
  write.table(out_dn, paste0("analysis/integration/summary/", fname_dn),
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              sep = ",")
}

# Immune sigs
for(i in 1:length(imm_n_gsl)){
  id <- names(imm_n_gsl)[i]
  
  out_gene <- data.frame(gene = imm_n_gsl[[i]]$gene)
  fname_gene <- paste0(id, "_sig_genes.csv")
  write.table(out_gene, paste0("analysis/integration/summary/", fname_gene),
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              sep = ",")
}