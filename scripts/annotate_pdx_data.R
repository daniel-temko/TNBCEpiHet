source("../utilities/utils.R")
source("../utilities/signatures.R")
source("../src/hexplot/hexplot.R")

setwd("../")

# Load data --------------------------------------------------------------------------------

load("R_Data/PDX_RNA_Data.RData")
load("R_Data/PDX_Heatmap_Metadata.RData")
load("R_Data/Gene_Sigs.RData")
tnbc_types <- c("basal", "mesenchymal", "luminal")
tnbc_type_names <- c(bas = "basal", mes = "mesenchymal", lum = "luminal", unassigned = "unassigned")

# Create gene signatures and assess in samples -------------------------------------------------------------------------

# Center data
pdx_rna_sig_data <- sweep(pdx_rna_rpkm_norm_uf, 1, rowMeans(pdx_rna_rpkm_norm_uf))

tnbc_gsl_pdx <- uniform_subset_gsl(tnbc_gsl, rownames(pdx_rna_sig_data))  

pdx_rna_sig_stats <- easy_signature_score(pdx_rna_sig_data, tnbc_gsl_pdx)
assigned_types <- score2type(pdx_rna_sig_stats)

# create the PDX cell line annotations
all(names(assigned_types) %in% pdx_cell_lines$sample)
pdx_line_annotation <- data.frame(row.names = pdx_cell_lines$sample, type = rep("not_available", nrow(pdx_cell_lines)))
pdx_line_annotation[names(assigned_types), "type"] <- tnbc_type_names[assigned_types]
pdx_line_annotation[names(assigned_types), "type2"] <- assigned_types

save(pdx_line_colours, pdx_line_annotation, file = "R_Data/PDX_Heatmap_Metadata_Annotated.RData")
