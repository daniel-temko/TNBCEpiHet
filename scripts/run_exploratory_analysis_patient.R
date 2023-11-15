
setwd("../")

source("scripts/utilities/exploratory_plots.R")
source("scripts/utilities/preprocessing.R")

# Load  data ----------------------------------------------------------------------------

load(file = "R_Data/Methyl_Data.RData")
load(file = "R_Data/Pat_Methyl_Data.RData")
load(file = "R_Data/Heatmap_Metadata.RData")

patient_colour <- "#800080"

# DNA Methylation data ----------------------------------------------------------------------------

pat_comb_methyl_metadata <- data.frame(row.names = c(rownames(methyl_metadata), rownames(pat_methyl_metadata)),
                                       type = c(methyl_metadata$tnbc.subtype, rep("patient", nrow(pat_methyl_metadata))))
pat_comb_methyl_ann_col <- pat_comb_methyl_metadata
pat_comb_methyl_ann_colours <- list(type = c(line_colours$type, patient = patient_colour))

pat_common_methyl_se <- intersect(rownames(methyl_se_agg_m), rownames(pat_methyl_cl_se_agg_m))
pat_comb_methyl_se_agg_m <- cbind(methyl_se_agg_m[pat_common_methyl_se, ], 
                                  pat_methyl_cl_se_agg_m[pat_common_methyl_se, ])

methyl_se_ranks <- variable_features(methyl_se_agg_m)
methyl_se_ranks <- methyl_se_ranks[methyl_se_ranks %in% pat_common_methyl_se]

GenerateExploratoryPlots(data = pat_comb_methyl_se_agg_m,
                         top_features = methyl_se_ranks,
                         n_pcs = rep(10, 6),
                         base_col = "#7C00FF",
                         base_dir = "analysis/patient_dna_methylation/exploratory/cl_superenhancer/",
                         annotation_col = pat_comb_methyl_ann_col,
                         annotation_colours = pat_comb_methyl_ann_colours)
