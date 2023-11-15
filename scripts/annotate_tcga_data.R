source("../utilities/utils.R")
source("../utilities/signatures.R")

setwd("../")

# Load data --------------------------------------------------------------------------------

load("R_Data/TCGA_RNA_Data.RData")
load("R_Data/TCGA_Heatmap_Metadata_Filtered.RData")
load("R_Data/Gene_Sigs.RData")
load("R_Data/PRRX1_Targets.RData")

tnbc_types <- c("bas", "mes", "lum")
tnbc_type_names <- c(bas = "basal", mes = "mesenchymal", lum = "luminal")

# Create gene signatures and assess in samples -------------------------------------------------------------------------

# Center data
tcga_tnbc_sig_data <- sweep(tcga_tnbc_rna_fpkm_norm_uf, 1, rowMeans(tcga_tnbc_rna_fpkm_norm_uf))

tnbc_gsl_tcga <- uniform_subset_gsl(tnbc_gsl, rownames(tcga_tnbc_sig_data))

tcga_sig_stats <- easy_signature_score(tcga_tnbc_sig_data, tnbc_gsl_tcga)
assigned_types <- score2type(tcga_sig_stats)
tcga_tnbc_rna_metadata$tnbc.type <- factor(assigned_types, levels = c(tnbc_types))

# Assign 'not_avaiable' to any sample without RNA data for type assignment (in this case there are none)
tcga_tnbc_sample_annotation$type <- "not_available"
tcga_tnbc_sample_annotation[names(assigned_types), "type"] <- tnbc_type_names[assigned_types]
tcga_tnbc_sample_annotation[names(assigned_types), "type2"] <- assigned_types

# PRRX1 expression and quantiles
#tcga_tnbc_prrx1_exp <- as.numeric(tcga_tnbc_rna_fpkm_norm_uf["PRRX1",])
tcga_tnbc_prrx1_exp <- as.numeric(tcga_tnbc_sig_data["PRRX1",])
tcga_tnbc_rna_metadata$prrx1 <- tcga_tnbc_prrx1_exp
tcga_tnbc_rna_metadata$prrx1.quantile <- get_quantiles(tcga_tnbc_prrx1_exp, 2)

# PRRX1 signature expression and quantiles
prrx1_gsl_tcga <- subset_gsl(gsl = prrx1_gsl, genes = rownames(tcga_tnbc_sig_data))

prrx1_sig_stats <- easy_signature_score(exp_data = tcga_tnbc_sig_data, gsl = prrx1_gsl_tcga)

tcga_tnbc_rna_metadata <- cbind(tcga_tnbc_rna_metadata, prrx1_sig_stats)

prrx1_tgt_quantiles <- lapply(as.data.frame(prrx1_sig_stats), function(x) get_quantiles(x, 2))
prrx1_tgt_quantiles <- as.data.frame(prrx1_tgt_quantiles)
colnames(prrx1_tgt_quantiles) <- paste0(colnames(prrx1_sig_stats), ".quantile")
tcga_tnbc_rna_metadata <- cbind(tcga_tnbc_rna_metadata, prrx1_tgt_quantiles)

# Assess immune signature scores in TCGA TNBC samples ----------------------------------------------------------

imm_n_gsl_tcga <- subset_gsl(gsl = imm_n_gsl, genes = rownames(tcga_tnbc_sig_data))

imm_n_sig_scores_tcga <- sapply(imm_n_gsl_tcga, function(x) {
  colMeans(tcga_tnbc_sig_data[x$gene, ])
})

imm_n_sig_quantiles_tcga <- lapply(as.data.frame(imm_n_sig_scores_tcga), function(x) get_quantiles(x, 2))
imm_n_sig_quantiles_tcga <- as.data.frame(imm_n_sig_quantiles_tcga)
rownames(imm_n_sig_quantiles_tcga) <- rownames(imm_n_sig_scores_tcga)
colnames(imm_n_sig_quantiles_tcga) <- paste0(colnames(imm_n_sig_quantiles_tcga), ".quantile")

imm_gsl_tcga <- subset_gsl(gsl = imm_gsl, genes = rownames(tcga_tnbc_sig_data))

imm_sig_scores_tcga <- sapply(imm_gsl_tcga, function(x){
  colMeans(tcga_tnbc_sig_data[x$gene,])
})

imm_sig_quantiles_tcga <- lapply(as.data.frame(imm_sig_scores_tcga), function(x) get_quantiles(x, 2))
imm_sig_quantiles_tcga <- as.data.frame(imm_sig_quantiles_tcga)
rownames(imm_sig_quantiles_tcga) <- rownames(imm_sig_scores_tcga)
colnames(imm_sig_quantiles_tcga) <- paste0(colnames(imm_sig_quantiles_tcga), ".quantile")

# Save results ------------------------------------------------------------------------------------------

save(tcga_rna_fpkm_norm_uf,
     tcga_rna_metadata,
     tcga_rna_patient_metadata,
     tcga_rna_uf_ann,
     tcga_tnbc_rna_fpkm_norm_uf,
     tcga_tnbc_rna_metadata,
     tcga_tnbc_rna_patient_metadata,
     imm_n_sig_scores_tcga,
     imm_n_sig_quantiles_tcga,
     imm_sig_scores_tcga,
     imm_sig_quantiles_tcga,
     file = "R_Data/TCGA_RNA_Annotated.RData")
save(tcga_tnbc_line_colours,
     tcga_sample_annotation,
     tcga_patient_data,
     tcga_tnbc_sample_annotation,
     tcga_tnbc_patient_data,
     file = "R_Data/TCGA_Heatmap_Metadata_Annotated.RData")
