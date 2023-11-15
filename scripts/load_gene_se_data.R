source("../utilities/utils.R")

setwd("../")

load("R_Data/Methyl_Data.RData")
load("R_Data/ChIP_Data.RData")

# Create gene-based H3K27ac ChIP-seq dataset
chip_se_ann$region.midpoint <- round((chip_se_ann$end.region + chip_se_ann$start.region) / 2)
chip_se_ann$dist.to.tss <- abs(chip_se_ann$region.midpoint - chip_se_ann$start.tss)
chip_se_ann$dist.to.tss[which(chip_se_ann$gene.symbol == ".")] <- -1

reg_genes <- unique(chip_se_ann$gene.symbol)
reg_genes <- reg_genes[which(reg_genes != ".")]

#' Create mapping from genes to top-regulating super-enhancer. Retain only
#' the super-enhancers that regulate a gene -- i.e. exclude gene.symbol == "." 
#' -- and take the closest super-enhancer in cases where a gene is nearest to
#' more than one super-enhancer
gene_se_map <- chip_se_ann[which(chip_se_ann$gene.symbol %in% reg_genes), ]
gene_se_map <- gene_se_map[order(gene_se_map$dist.to.tss), ]
gene_se_map <- gene_se_map[match(reg_genes, gene_se_map$gene.symbol),]
gene_se_map <- rownames_from_col(bring_to_front(gene_se_map, "gene.symbol"))

# Create chip_gene_se_data
keep <- gene_se_map$dim.id %in% rownames(chip_se_all_rpkm_norm_uf)
chip_gene_se_ann <- gene_se_map[keep, ]
chip_gene_se_rpkm_norm_uf <- chip_se_all_rpkm_norm_uf[as.character(chip_gene_se_ann$dim.id), ]
rownames(chip_gene_se_rpkm_norm_uf) <- rownames(chip_gene_se_ann)

# Create methyl_gene_se_int_agg_m
keep <- gene_se_map$dim.id %in% rownames(methyl_se_int_agg_m)
methyl_gene_se_int_agg_ann <- gene_se_map[keep, ]
methyl_gene_se_int_agg_m <- methyl_se_int_agg_m[as.character(methyl_gene_se_int_agg_ann$dim.id), ]
rownames(methyl_gene_se_int_agg_m) <- rownames(methyl_gene_se_int_agg_ann)

# Create methyl_gene_se_agg_m
keep <- gene_se_map$dim.id %in% rownames(methyl_se_agg_m)
methyl_gene_se_agg_ann <- gene_se_map[keep, ]
methyl_gene_se_agg_m <- methyl_se_agg_m[as.character(methyl_gene_se_agg_ann$dim.id), ]
rownames(methyl_gene_se_agg_m) <- rownames(methyl_gene_se_agg_ann)

save(chip_gene_se_rpkm_norm_uf,
     chip_gene_se_ann,
     methyl_gene_se_agg_m,
     methyl_gene_se_agg_ann,
     methyl_gene_se_int_agg_m, 
     methyl_gene_se_int_agg_ann,
     file = "R_Data/Gene_SE_Data.RData")