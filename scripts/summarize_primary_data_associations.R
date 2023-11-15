common_sig_list <- list()
for(i in 1:8){
  fid <- paste0("factor", i)
  df1 <- read.table(paste0("analysis/human_primary_data/metabric/mofa_factors/immune_signatures_and_mofa_factors/", fid, "/lists/Factor", i, "_correlation_test_p_values.csv"), sep = ",", header = TRUE)
  df2 <- read.table(paste0("analysis/human_primary_data/tcga/mofa_factors/immune_signatures_and_mofa_factors/", fid, "/lists/Factor", i, "_correlation_test_p_values.csv"), sep = ",", header = TRUE)  
  sig_up1 <- subset(df1, p.adj < 0.05 & rho > 0)
  sig_up2 <- subset(df2, p.adj < 0.05 & rho > 0)
  sig_dn1 <- subset(df1, p.adj < 0.05 & rho < 0)
  sig_dn2 <- subset(df2, p.adj < 0.05 & rho < 0)
  message("F", i)
  message("MBC sig: ", nrow(sig1))
  message("TCGA sig: ", nrow(sig2))
  common_up_sig <- intersect(sig_up1$sig.id, sig_up2$sig.id)
  common_dn_sig <- intersect(sig_dn1$sig.id, sig_dn2$sig.id)
  message("Both up:", length(common_up_sig), " Both dn:", length(common_dn_sig))
  common_sig_list[[fid]] <- list(up = common_up_sig, dn = common_dn_sig)
}

