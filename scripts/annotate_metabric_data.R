options(warn = 2)
source("../utilities/utils.R")
source("../utilities/mofa2.R")
source("../utilities/signatures.R")

setwd("../")

# Load data --------------------------------------------------------------------------------

load("R_Data/Metabric_Data.RData")
load("R_Data/Gene_Sigs.RData")
load("R_Data/PRRX1_Targets.RData")
load("R_Data/MOFA2_TopN_Heatmap_Feature_Filtering_Seed20200225_RNA_Weights.RData")

tnbc_types <- c("bas", "mes", "lum")

# Analysis ---------------------------------------------------------------------------

# Center
mbc_tnbc_sig_data <- sweep(mbc_tnbc_rna_med, 1, rowMeans(mbc_tnbc_rna_med))

# TNBC type signatures
tnbc_gsl_mbc <- uniform_subset_gsl(tnbc_gsl, rownames(mbc_tnbc_sig_data))

mbc_sig_stats <- easy_signature_score(mbc_tnbc_sig_data, tnbc_gsl_mbc)
assigned_types <- score2type(mbc_sig_stats)
mbc_tnbc_metadata$tnbc.type <- factor(assigned_types, levels = tnbc_types)

# PRRX1 expression and quantiles
mbc_tnbc_prrx1_exp <- as.numeric(mbc_tnbc_sig_data["PRRX1",])
mbc_tnbc_metadata$prrx1 <- mbc_tnbc_prrx1_exp
mbc_tnbc_metadata$prrx1.quantile <- get_quantiles(mbc_tnbc_prrx1_exp, 2)

# PRRX1 signature expression and quantiles
prrx1_gsl_mbc <- subset_gsl(prrx1_gsl, rownames(mbc_tnbc_sig_data))

prrx1_sig_stats <- easy_signature_score(exp_data = mbc_tnbc_sig_data, gsl = prrx1_gsl_mbc)

mbc_tnbc_metadata <- cbind(mbc_tnbc_metadata, prrx1_sig_stats)

prrx1_tgt_quantiles <- lapply(as.data.frame(prrx1_sig_stats), function(x) get_quantiles(x, 2))
prrx1_tgt_quantiles <- as.data.frame(prrx1_tgt_quantiles)
colnames(prrx1_tgt_quantiles) <- paste0(colnames(prrx1_tgt_quantiles), ".quantile")
mbc_tnbc_metadata <- cbind(mbc_tnbc_metadata, prrx1_tgt_quantiles)

# Assess immune signature scores in Metabric TNBC samples 
imm_n_gsl_mbc <- subset_gsl(imm_n_gsl, rownames(mbc_tnbc_sig_data))

imm_n_sig_scores_mbc <- sapply(imm_n_gsl_mbc, function(x){
  colMeans(mbc_tnbc_sig_data[x$gene,])
})

imm_n_sig_quantiles_mbc <- lapply(as.data.frame(imm_n_sig_scores_mbc), function(x) get_quantiles(x, 2))
imm_n_sig_quantiles_mbc <- as.data.frame(imm_n_sig_quantiles_mbc)
rownames(imm_n_sig_quantiles_mbc) <- rownames(imm_n_sig_scores_mbc)
colnames(imm_n_sig_quantiles_mbc) <- paste0(colnames(imm_n_sig_quantiles_mbc), ".quantile")

imm_gsl_mbc <- subset_gsl(imm_gsl, rownames(mbc_tnbc_sig_data))

imm_sig_scores_mbc <- sapply(imm_gsl_mbc, function(x){
  colMeans(mbc_tnbc_sig_data[x$gene,])
})

imm_sig_quantiles_mbc <- lapply(as.data.frame(imm_sig_scores_mbc), function(x) get_quantiles(x, 2))
imm_sig_quantiles_mbc <- as.data.frame(imm_sig_quantiles_mbc)
rownames(imm_sig_quantiles_mbc) <- rownames(imm_sig_scores_mbc)
colnames(imm_sig_quantiles_mbc) <- paste0(colnames(imm_sig_quantiles_mbc), ".quantile")

# Save results ------------------------------------------------------------------------------------------

save(mbc_rna_med,
     mbc_metadata,
     mbc_patient_metadata,
     mbc_tnbc_rna_med,
     mbc_tnbc_metadata,
     mbc_tnbc_patient_metadata,
     imm_n_sig_scores_mbc,
     imm_n_sig_quantiles_mbc,
     imm_sig_scores_mbc,
     imm_sig_quantiles_mbc,
     file = "R_Data/Metabric_Annotated.RData")
