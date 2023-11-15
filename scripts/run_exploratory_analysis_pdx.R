
setwd("../")

source("scripts/utilities/exploratory_plots.R")
source("scripts/utilities/preprocessing.R")

# Load  data ----------------------------------------------------------------------------

load(file = "R_Data/RNA_Data.RData")
load(file = "R_Data/PDX_RNA_Data.RData")
load(file = "R_Data/ChIP_Data.RData")
load(file = "R_Data/PDX_ChIP_Data.RData")
load(file = "R_Data/Methyl_Data.RData")
load(file = "R_Data/PDX_Methyl_Data.RData")
load(file = "R_Data/Gene_Sigs.RData")
load(file = "R_Data/Heatmap_Metadata.RData")

pdx_colour <- "#800080"

# Analyse RNA-seq data ----------------------------------------------------------------------------

comb_rna_metadata <- data.frame(row.names = c(rownames(rna_metadata), rownames(pdx_rna_metadata)),
                                type = c(rna_metadata$tnbc.subtype, rep("pdx", nrow(pdx_rna_metadata))))
comb_rna_ann_col <- comb_rna_metadata
comb_rna_ann_colours <- list(type = c(line_colours$type, pdx = pdx_colour))

all(sort(rownames(rna_rpkm_uf)) == sort(rownames(pdx_rna_rpkm_uf)))
comb_rna_rpkm_uf <- cbind(rna_rpkm_uf,
                          pdx_rna_rpkm_uf[rownames(rna_rpkm_uf), ])

all(sort(rownames(rna_rpkm_norm_uf)) == sort(rownames(pdx_rna_rpkm_norm_uf)))
comb_rna_rpkm_norm_uf <- cbind(rna_rpkm_norm_uf,
                               pdx_rna_rpkm_norm_uf[rownames(rna_rpkm_norm_uf), ])

# Most variable genes across the data
keep <- apply(comb_rna_rpkm_uf, 1, function(x) length(which(x > 1))) > 1
comb_rna_rpkm <- comb_rna_rpkm_uf[keep, ]
all(apply(comb_rna_rpkm, 1, function(x) var(x) > 0))

comb_rna_rpkm_norm <- log2(comb_rna_rpkm + 1)

comb_rna_ranks <- variable_features(comb_rna_rpkm_norm)

GenerateExploratoryPlots(data = comb_rna_rpkm_norm,
                         top_features = comb_rna_ranks,
                         n_pcs = rep(10, 6),
                         base_col = "#FF00FF",
                         base_dir = "analysis/pdx_rna_seq/exploratory/",
                         annotation_col = comb_rna_ann_col,
                         annotation_colours = comb_rna_ann_colours)

# Cell line genes
rna_ranks <- variable_features(rna_rpkm_norm)

GenerateExploratoryPlots(data = comb_rna_rpkm_norm_uf[rownames(rna_rpkm_norm), ],
                         top_features = rna_ranks,
                         n_pcs = rep(10, 6),
                         base_col = "#FF00FF",
                         base_dir = "analysis/pdx_rna_seq/exploratory_cl_genes/",
                         annotation_col = comb_rna_ann_col,
                         annotation_colours = comb_rna_ann_colours)

# Cell line TNBC genes
bas_genes <- unique(unlist(tnbc_gsl$bas))
lum_genes <- unique(unlist(tnbc_gsl$lum))
mes_genes <- unique(unlist(tnbc_gsl$mes))
tnbc_genes <- unique(c(bas_genes, lum_genes, mes_genes))
all(tnbc_genes %in% rownames(rna_rpkm_norm_uf))
rna_rpkm_norm_tnbc_genes <- rna_rpkm_norm_uf[tnbc_genes, ]
rna_tnbc_ranks <- variable_features(rna_rpkm_norm_tnbc_genes)

GenerateExploratoryPlots(data = comb_rna_rpkm_norm_uf[tnbc_genes, ],
                         top_features = rna_tnbc_ranks,
                         n_pcs = rep(10, 6),
                         base_col = "#FF00FF",
                         base_dir = "analysis/pdx_rna_seq/exploratory_cl_tnbc_genes/",
                         annotation_col = comb_rna_ann_col,
                         annotation_colours = comb_rna_ann_colours)

# Analyse ChIP-seq data ----------------------------------------------------------------------------

comb_chip_metadata <- data.frame(row.names = c(rownames(chip_metadata), rownames(pdx_chip_metadata)),
                                 type = c(chip_metadata$tnbc.subtype, rep("pdx", nrow(pdx_chip_metadata))))
comb_chip_ann_col <- comb_chip_metadata
comb_chip_ann_colours <- list(type = c(line_colours$type, pdx = pdx_colour))

all(sort(rownames(chip_se_all_rpkm_norm_uf)) == sort(rownames(pdx_chip_cl_se_all_rpkm_norm_uf)))
comb_chip_cl_se_all_rpkm_norm_uf <- cbind(chip_se_all_rpkm_norm_uf,
                                          pdx_chip_cl_se_all_rpkm_norm_uf[rownames(chip_se_all_rpkm_norm_uf), ])

chip_ranks <- variable_features(chip_se_all_rpkm_norm)

GenerateExploratoryPlots(data = comb_chip_cl_se_all_rpkm_norm_uf[rownames(chip_se_all_rpkm_norm), ],
                         top_features = chip_ranks,
                         n_pcs = rep(10, 6),
                         base_col = "#FF8000",
                         base_dir = "analysis/pdx_chip_seq/exploratory_cl_superenhancers/",
                         annotation_col = comb_chip_ann_col,
                         annotation_colours = comb_chip_ann_colours)

# DNA Methylation data ----------------------------------------------------------------------------

comb_methyl_metadata <- data.frame(row.names = c(rownames(methyl_metadata), rownames(pdx_methyl_metadata)),
                                   type = c(methyl_metadata$tnbc.subtype, rep("pdx", nrow(pdx_methyl_metadata))))
comb_methyl_ann_col <- comb_methyl_metadata
comb_methyl_ann_colours <- list(type = c(line_colours$type, pdx = pdx_colour))

common_methyl_se <- intersect(rownames(methyl_se_agg_m), rownames(pdx_methyl_cl_se_agg_m))
comb_methyl_se_agg_m <- cbind(methyl_se_agg_m[common_methyl_se, ], pdx_methyl_cl_se_agg_m[common_methyl_se, ])

methyl_se_ranks <- variable_features(methyl_se_agg_m)
methyl_se_ranks <- methyl_se_ranks[methyl_se_ranks %in% common_methyl_se]

GenerateExploratoryPlots(data = comb_methyl_se_agg_m,
                         top_features = methyl_se_ranks,
                         n_pcs = rep(10, 6),
                         base_col = "#7C00FF",
                         base_dir = "analysis/pdx_dna_methylation/exploratory/cl_superenhancer/",
                         annotation_col = comb_methyl_ann_col,
                         annotation_colours = comb_methyl_ann_colours)
