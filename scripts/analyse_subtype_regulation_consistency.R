setwd("../")

################################################################################
# load library
################################################################################

library(pheatmap)
library(reshape)
library(ggplot2)
library(ggpubr)
library(FSA)
source("../utilities/utils.R")

################################################################################
# local functions
################################################################################

up_tfs <- function(df){
  df[which((df$log2FoldChange > 0) & (df$dbd.tf)),]
}

local_plot_boxplot <- function(df, type, ymin, ymax, y_lab, pal=line_colours$type){
    gg <- melt(t(df))
    colnames(gg) <- c("Cell.Line", "TF", "Expression")
    gg$type <- type
    
    ggplot(gg, aes(x=TF, y=Expression, fill=type)) +
        my_theme() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ylim(c(ymin, ymax)) +
        ylab(y_lab) +
        geom_boxplot(show.legend=FALSE, outlier.shape=NA) + # to prevent double plotting outliers
        geom_jitter(show.legend=FALSE, width=0.2) +
        geom_hline(yintercept=0, linetype="dashed", color="grey", size=1.5) +
        scale_fill_manual(values = pal)
}

local_plot_kw_boxplot <- function(gg, plot_col, y_lab, p_crit=0.05, pal=line_colours$type){
    gg$Value <- gg[[plot_col]]
    kt <- kruskal.test(Value ~ Type, data = gg)
    op <- round(kt$p.value, 3)
    p <- ggboxplot(gg, x="Type", y="Value", add=c("jitter"), fill="Type", 
                   ylab=y_lab, palette=pal)
    y_p <- max(gg$Value) * 1.15
    if(op < p_crit){
        dt <- dunnTest(Value ~ Type, data=gg, method="holm")
        g1 <- sapply(dt$res$Comparison, function(x) strsplit(x, " - ")[[1]][1])
        g2 <- sapply(dt$res$Comparison, function(x) strsplit(x, " - ")[[1]][2])
        st <- data.frame(group1 = g1, group2 = g2, p.adj = round(dt$res$P.adj, 3))
        y_max <- max(gg$Value) * 1.5
        p <- p + 
            stat_compare_means(label.y = y_max) +
            stat_pvalue_manual(st, label="p.adj", y.position=y_p, step.increase = 0.1)
    } else {
        p <- p +
            stat_compare_means(label.y = y_p)
    }
    return(p)
}

################################################################################
# load data
################################################################################

load("R_Data/RNA_Data.RData")
load("R_Data/ChIP_Data.RData")
load("R_Data/Heatmap_Metadata.RData")

bas_de <- read_csv("analysis/rna_seq/differential_genes/basal_de_lfc0_genes.csv")
mes_de <- read_csv("analysis/rna_seq/differential_genes/mesenchymal_de_lfc0_genes.csv")
lum_de <- read_csv("analysis/rna_seq/differential_genes/luminal_de_lfc0_genes.csv")

bas_dse <- read_csv("analysis/chip_seq/differential_regions/superenhancers_basal_diff_lfc0.csv")
mes_dse <- read_csv("analysis/chip_seq/differential_regions/superenhancers_mesenchymal_diff_lfc0.csv")
lum_dse <- read_csv("analysis/chip_seq/differential_regions/superenhancers_luminal_diff_lfc0.csv")

bas_tfs <- intersect(up_tfs(bas_de)$gene.symbol, up_tfs(bas_dse)$closest.gene)
mes_tfs <- intersect(up_tfs(mes_de)$gene.symbol, up_tfs(mes_dse)$closest.gene)
lum_tfs <- intersect(up_tfs(lum_de)$gene.symbol, up_tfs(lum_dse)$closest.gene)

bas_ses <- subset(up_tfs(bas_dse), closest.gene %in% bas_tfs)$region
mes_ses <- subset(up_tfs(mes_dse), closest.gene %in% mes_tfs)$region
lum_ses <- subset(up_tfs(lum_dse), closest.gene %in% lum_tfs)$region

################################################################################
# transform data
################################################################################

stopifnot(all(c(bas_ses, mes_ses, lum_ses) %in% rownames(chip_se_all_rpkm_norm_uf)))
stopifnot(all(c(bas_tfs, mes_tfs, lum_tfs) %in% rownames(rna_rpkm_norm_uf)))

# Form RNA and ChIP matrices
rna_data <- rna_rpkm_norm_uf[c(bas_tfs, mes_tfs, lum_tfs), ]
rna_types <- line_annotation[colnames(rna_data), "type"]
rna_data <- rna_data[,order(rna_types)]
rna_types <- sort(rna_types)
rna_data_centered <- sweep(rna_data, 1, rowMeans(rna_data))

chip_data <- chip_se_all_rpkm_norm_uf[c(bas_ses, mes_ses, lum_ses), ]
chip_types <- line_annotation[colnames(chip_data), "type"]
chip_data <- chip_data[,order(chip_types)]
chip_types <- sort(chip_types)
chip_data_centered <- sweep(chip_data, 1, rowMeans(chip_data))

################################################################################
# ChIP heatmap
################################################################################

# ordered version
pp <- chip_data_centered

bas_se_ids <- rename_duplicates(chip_se_ann[bas_ses, "gene.symbol"])
mes_se_ids <- rename_duplicates(chip_se_ann[mes_ses, "gene.symbol"])
lum_se_ids <- rename_duplicates(chip_se_ann[lum_ses, "gene.symbol"])

rownames(pp) <- c(bas_se_ids, mes_se_ids, lum_se_ids)
pp <- pp[c(sort(bas_se_ids), sort(mes_se_ids), sort(lum_se_ids)),]

#chip_gene_ids <- chip_se_ann[c(bas_ses, mes_ses, lum_ses), "gene.symbol"]
#rownames(pp) <- paste0(chip_gene_ids, " ", rownames(pp))
#rownames(pp) <- rename_duplicates(chip_gene_ids)

pdf("analysis/integration/subtype_regulation/de_dse_tfs_chip_lfc.pdf")
  p <- pheatmap(pp, 
           annotation_col = line_annotation, 
           annotation_colors = line_colours, 
           cluster_rows = FALSE,
           cluster_cols = FALSE)
  print(p)
dev.off()

# clustered version
pp <- chip_data_centered
rownames(pp) <- rename_duplicates(chip_se_ann[rownames(pp), "gene.symbol"])

col_dist <- dist(t(pp), method = "euclidean")
col_hclust <- hclust(col_dist, method = "ward.D2")

row_dist <- dist(pp, method = "euclidean")
row_hclust <- hclust(row_dist, method = "ward.D2")

se_annotation <- data.frame(row.names = rownames(pp),
                            type = c(rep("basal", length(bas_ses)),
                                     rep("mesenchymal", length(mes_ses)),
                                     rep("luminal", length(lum_ses))))

pdf("analysis/integration/subtype_regulation/de_dse_tfs_chip_lfc_clust.pdf")
p <- pheatmap(pp, 
              annotation_col = line_annotation, 
              annotation_row = se_annotation,
              annotation_colors = line_colours, 
              cluster_rows = row_hclust,
              cluster_cols = col_hclust)
print(p)
dev.off()

################################################################################
# RNA variability
################################################################################

# Analyse consistency of expression
rna_bas <- rna_data_centered[bas_tfs, which(rna_types == "basal")]
rna_lum <- rna_data_centered[lum_tfs, which(rna_types == "luminal")]
rna_mes <- rna_data_centered[mes_tfs, which(rna_types == "mesenchymal")]

pdf("analysis/integration/subtype_regulation/tf_expression_dist_bas.pdf", 9, 9)
local_plot_boxplot(rna_bas, "basal", ymin=-4, ymax=6.5, y_lab="Expression")
dev.off()

pdf("analysis/integration/subtype_regulation/tf_expression_dist_lum.pdf", 9, 9)
local_plot_boxplot(rna_lum, "luminal", ymin=-4, ymax=6.5, y_lab="Expression")
dev.off()

pdf("analysis/integration/subtype_regulation/tf_expression_dist_mes.pdf", 9, 9)
local_plot_boxplot(rna_mes, "mesenchymal", ymin=-4, ymax=6.5, y_lab="Expression")
dev.off()

################################################################################
# H3K27ac variability
################################################################################

chip_bas <- chip_data_centered[bas_ses, which(chip_types == "basal")]
chip_lum <- chip_data_centered[lum_ses, which(chip_types == "luminal")]
chip_mes <- chip_data_centered[mes_ses, which(chip_types == "mesenchymal")]

rownames(chip_bas) <- rename_duplicates(chip_se_ann[bas_ses, "gene.symbol"])
rownames(chip_lum) <- rename_duplicates(chip_se_ann[lum_ses, "gene.symbol"])
rownames(chip_mes) <- rename_duplicates(chip_se_ann[mes_ses, "gene.symbol"])

pdf("analysis/integration/subtype_regulation/tf_h3k27ac_dist_bas.pdf", 9, 9)
local_plot_boxplot(chip_bas, "basal", ymin=-3, ymax=5, y_lab="H3K27ac Expression")
dev.off()

pdf("analysis/integration/subtype_regulation/tf_h3k27ac_dist_lum.pdf", 9, 9)
local_plot_boxplot(chip_lum, "luminal", ymin=-3, ymax=5, y_lab="H3K27ac Expression")
dev.off()

pdf("analysis/integration/subtype_regulation/tf_h3k27ac_dist_mes.pdf", 9, 9)
local_plot_boxplot(chip_mes, "mesenchymal", ymin=-3, ymax=5, y_lab="H3K27ac Expression")
dev.off()

################################################################################
# Compare variances
################################################################################

# get uncentered versions of the subsetted matrices
rna_bas <- rna_data[bas_tfs, which(rna_types == "basal")]
rna_lum <- rna_data[lum_tfs, which(rna_types == "luminal")]
rna_mes <- rna_data[mes_tfs, which(rna_types == "mesenchymal")]

rna_var_bas <- apply(rna_bas, 1, var)
rna_var_lum <- apply(rna_lum, 1, var)
rna_var_mes <- apply(rna_mes, 1, var)

rna_mean_bas <- rowMeans(rna_bas)
rna_mean_lum <- rowMeans(rna_lum)
rna_mean_mes <- rowMeans(rna_mes)

gg <- data.frame(TF = c(bas_tfs, lum_tfs, mes_tfs),
                 Var = c(rna_var_bas, rna_var_lum, rna_var_mes),
                 Mean = c(rna_mean_bas, rna_mean_lum, rna_mean_mes),
                 Type = c(rep("basal", length(bas_tfs)), 
                          rep("luminal", length(lum_tfs)), 
                          rep("mesenchymal", length(mes_tfs))))
gg$Type <- factor(gg$Type)

pdf("analysis/integration/subtype_regulation/rna_var_summary.pdf", 5, 7)
local_plot_kw_boxplot(gg, "Var", "Expression Variance")
dev.off()

pdf("analysis/integration/subtype_regulation/rna_mean_summary.pdf", 5, 7)
local_plot_kw_boxplot(gg, "Mean", "Mean Expression")
dev.off()

# get uncentered versions of the subsetted matrices
chip_bas <- chip_data[bas_ses, which(chip_types == "basal")]
chip_lum <- chip_data[lum_ses, which(chip_types == "luminal")]
chip_mes <- chip_data[mes_ses, which(chip_types == "mesenchymal")]

rownames(chip_bas) <- rename_duplicates(chip_se_ann[bas_ses, "gene.symbol"])
rownames(chip_lum) <- rename_duplicates(chip_se_ann[lum_ses, "gene.symbol"])
rownames(chip_mes) <- rename_duplicates(chip_se_ann[mes_ses, "gene.symbol"])

chip_var_bas <- apply(chip_bas, 1, var)
chip_var_lum <- apply(chip_lum, 1, var)
chip_var_mes <- apply(chip_mes, 1, var)

chip_mean_bas <- rowMeans(chip_bas)
chip_mean_lum <- rowMeans(chip_lum)
chip_mean_mes <- rowMeans(chip_mes)

gg <- data.frame(SE = c(bas_ses, lum_ses, mes_ses),
                 Var = c(chip_var_bas, chip_var_lum, chip_var_mes),
                 Mean = c(chip_mean_bas, chip_mean_lum, chip_mean_mes),
                 Type = c(rep("basal", length(bas_ses)), 
                          rep("luminal", length(lum_ses)), 
                          rep("mesenchymal", length(mes_ses))))

gg$Type <- factor(gg$Type)

pdf("analysis/integration/subtype_regulation/h3k27ac_var_summary.pdf", 5, 7)
local_plot_kw_boxplot(gg, "Var", "H3K27ac Expression Variance")
dev.off()

pdf("analysis/integration/subtype_regulation/h3k27ac_mean_summary.pdf", 5, 7)
local_plot_kw_boxplot(gg, "Mean", "Mean H3K27ac Expression")
dev.off()
