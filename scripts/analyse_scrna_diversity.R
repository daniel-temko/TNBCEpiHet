#####
# analyze each sample using Negative Binomial distribution
#####

setwd("../")

################################################################################
# load library
################################################################################
source("scripts/utilities/utils.R")
library(reshape)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(FSA)
library(ggpubr)
library(MASS)
library(lme4)
library(splines)
library(car)
library(forestplot)

################################################################################
# load data
################################################################################
obj_phe <- load("R_Data/SC_Metadata.RData")
obj_met <- load("R_Data/Heatmap_Metadata.RData")
obj_dat <- load("R_Data/SC_Filt.RData")
obj_sig <- load("R_Data/Gene_Sigs.RData")
obj_deg <- load("R_Data/Subtype_DE_Genes.RData")
obj_rna <- load("R_Data/RNA_Data.RData")

hkg <- c("ACTB", "GAPDH", "C1orf43", "CHMP2A", "GPI", "PSMB2", "PSMB4", "RAB7A",
         "REEP5", "SNRPD3", "VCP", "VPS29")

min_exp_val <- 1
scale_fac <- 10000

cl_samp_num <- rowSums(table(scrna_sample_table$cell.line, scrna_sample_table$chemistry))
rep_cl <- names(cl_samp_num)[which(cl_samp_num > 1)]
test_exc_samples <- paste0(rep_cl, "_B2")

################################################################################
# transform data
################################################################################

## filter single-cell data
exp_num <- rowSums(sc[["RNA"]]@counts > 0)
keep <- exp_num >= 10
sc <- sc[keep, ]

################################################################################
# local functions
################################################################################

est_scv <- function(exp_count_mat, scale_factor = scale_fac){
  size_factors <- colSums(exp_count_mat) / scale_factor
  exp_norm_mat <- sweep(exp_count_mat, 2, size_factors, "/")
  gene_cs_vars <- apply(exp_norm_mat, 1, var)
  gene_cs_means <- rowMeans(exp_norm_mat)
  bias_factors <- (gene_cs_means / ncol(exp_count_mat)) * sum(1 / size_factors)
  est_raw_var <- pmax(gene_cs_vars - bias_factors, 0)
  #est_raw_var <- gene_cs_vars - bias_factors
  est_scv <- est_raw_var / gene_cs_means^2
  list(cells = data.frame(size_factors = size_factors),
       genes = data.frame(cs_mean = gene_cs_means, 
                          cs_var = gene_cs_vars, 
                          est_raw_var = est_raw_var,
                          est_scv = est_scv,
                          est_theta = 1 / est_scv))
}

est_nb_reg <- function(exp_count_mat, gene_names, scale_factor = scale_fac){
  size_factors <- colSums(exp_count_mat) / scale_factor
  exp_norm_mat <- sweep(exp_count_mat, 2, size_factors, "/")
  gene_cs_means <- rowMeans(exp_norm_mat)
  ## subset to required genes only
  gene_cs_means <- gene_cs_means[gene_names]
  exp_count_mat <- exp_count_mat[gene_names, ]
  theta_vals <- sapply(1:nrow(exp_count_mat), function(x){
    if(x %% 50 == 0) message(x, " of ", nrow(exp_count_mat))
    count_vals <- exp_count_mat[x, ]
    stopifnot(sum(count_vals) > 0)
    tryCatch(
      warning = function(cnd) NA,
      {
        offs <- log(gene_cs_means[x] * size_factors)
        res <- glm.nb(count_vals ~ offset(offs) + 0)
        res$theta  
      }
    )
  })
  list(cells = data.frame(size_factors = size_factors),
       genes = data.frame(row.names = gene_names,
                          cs_mean = gene_cs_means,
                          est_theta = theta_vals))
}

local_boxplot <- function(f_xx, y_col, title_str, y_str, pal = line_colours$type, cap_val = Inf, 
                          p_cut = 0.05, asinh_y = FALSE, asinh_cofactor = 5, log_y = FALSE,
                          log_pseudocount = 0){
  stopifnot(is.factor(f_xx$tnbc_type) & 
              identical(levels(f_xx$tnbc_type), c("basal", "luminal", "mesenchymal")))
  if(asinh_y){
    f_xx$y <- asinh(f_xx[[y_col]] * asinh_cofactor)
    y_lab <- paste0("asinh(", y_str, " * ", asinh_cofactor, ")")
  } else if(log_y) {
    f_xx$y <- log(f_xx[[y_col]] + log_pseudocount)
    y_lab <- paste0("log(", y_str, " + ", log_pseudocount, ")")
  } else {
    f_xx$y <- f_xx[[y_col]]  
    y_lab <- y_str
  }
  if(!cap_val == Inf){
    f_xx$y <- pmin(f_xx$y, cap_val)  
  }
  kw_res <- kruskal.test(y ~ tnbc_type, f_xx)
  dunn_res <- dunnTest(y ~ tnbc_type, f_xx)
  y_bounds <- c(min(f_xx$y), max(f_xx$y))
  y_scale <- (y_bounds[2] - y_bounds[1])
  gg_s <- format_n(f_xx, y_bounds[2] + y_scale/10)
  p <- ggplot(f_xx) +
    geom_boxplot(aes(x = tnbc_type, y = y, fill = tnbc_type), outlier.shape = NA) +
    geom_jitter(aes(x = tnbc_type, y = y), color = "black", width = 0.1) +
    scale_fill_manual(values = pal) +
    ggtitle(title_str) +
    ylab(y_lab) +
    my_theme() +
    geom_text(data = gg_s,
              show.legend = FALSE,
              aes(x = x, y = y, label = n))
  if(kw_res$p.value < p_cut){
    gg_p <- format_dunn_p(dunn_res, y_bounds[2] + ((2:4) * (y_scale / 10)))
    p <- p + stat_pvalue_manual(gg_p, label = "p.adj")
  }
  print(p)
  return(list(kw = kw_res, dunn = dunn_res))
}

local_facet_boxplot <- function(f_xx, y_col, title_str, y_str, pal = line_colours$type, cap_val = Inf, 
                                p_cut = 0.05, asinh_y = FALSE, asinh_cofactor = 5, log_y = FALSE,
                                log_pseudocount = 0){
  stopifnot(is.factor(f_xx$tnbc_type) & 
              identical(levels(f_xx$tnbc_type), c("basal", "luminal", "mesenchymal")))
  if(asinh_y){
    f_xx$y <- asinh(f_xx[[y_col]] * asinh_cofactor)
    y_lab <- paste0("asinh(", y_str, " * ", asinh_cofactor, ")")
  } else if(log_y) {
    f_xx$y <- log(f_xx[[y_col]] + log_pseudocount)
    y_lab <- paste0("log(", y_str, " + ", log_pseudocount, ")")
  } else {
    f_xx$y <- f_xx[[y_col]]  
    y_lab <- y_str
  }
  if(!cap_val == Inf){
    f_xx$y <- pmin(f_xx$y, cap_val)  
  }
  y_bounds <- c(min(f_xx$y, na.rm = TRUE), max(f_xx$y, na.rm = TRUE))
  y_scale <- (y_bounds[2] - y_bounds[1])
  kw_list <- lapply(levels(f_xx$mean_bin), function(x){
    f_xx_sub <- subset(f_xx, mean_bin == x)
    if(length(which(is.na(f_xx_sub$y))) > 0){
      NA
    } else {
      kruskal.test(y ~ tnbc_type, f_xx_sub)
    }
  })
  dunn_list <- lapply(levels(f_xx$mean_bin), function(x){
    f_xx_sub <- subset(f_xx, mean_bin == x)
    if(length(which(is.na(f_xx_sub$y))) > 0){
      NA
    } else{
      dunn_res <- dunnTest(y ~ tnbc_type, f_xx_sub)
      gg_p <- cbind(mean_bin = x, format_dunn_p(dunn_res, max(f_xx_sub$y) + ((1:3) * (y_scale / 10))))
    }
  })
  p_ids <- which(sapply(kw_list, function(x){
    if(!class(x) == "htest"){
      FALSE
    } else{
      x$p.value < p_cut
    }
  }))
  p <- ggplot(f_xx) +
    geom_boxplot(aes(x = tnbc_type, y = y, fill = tnbc_type)) +
    scale_fill_manual(values = pal) +
    ggtitle(title_str) +
    ylab(y_lab) +
    my_theme() + 
    facet_grid(. ~ mean_bin)  
  if(length(p_ids) > 0){
    gg_p <- do.call(rbind, dunn_list[p_ids])
    p <- p + stat_pvalue_manual(gg_p, label = "p.adj")
  } 
  print(p)
  return(list(kw = kw_list, dunn = dunn_list))
}

format_n <- function(f_xx, y_val){
  gg_s <- data.frame(tnbc_type = levels(f_xx$tnbc_type))
  gg_s$x <- 1:nrow(gg_s)
  gg_s$n <- paste0("N=", sapply(gg_s$tnbc_type, function(x) length(which(f_xx$tnbc_type == x))))
  gg_s$y <- y_val
  gg_s
}

format_dunn_p <- function(dunn_res, y_vals, thresh = 0.05){
  sig_stars <- sapply(dunn_res$res$P.adj, function(x) {
    if(x < thresh) {
      "*"
    } else {
      ""
    }
  })
  gg_p <- data.frame(group1 = sapply(dunn_res$res$Comparison, function(x) strsplit(x, " - ")[[1]][1]),
                     group2 = sapply(dunn_res$res$Comparison, function(x) strsplit(x, " - ")[[1]][2]),
                     p.adj = paste0("P=", format(dunn_res$res$P.adj, scientific = TRUE, digits = 2), sig_stars),
                     y.position = y_vals)
}

plot_nb_fit <- function(exp_count_mat, gene_name, title_str, scale_factor = scale_fac){
  size_factors <- colSums(exp_count_mat) / scale_factor
  count_vals <- exp_count_mat[gene_name, ]
  cs_mean <- mean(count_vals / size_factors)
  offs <- log(cs_mean * size_factors)
  theta_est <- glm.nb(count_vals ~ offset(offs) + 0)$theta
  obs <- table(count_vals)
  density_list <- lapply(offs, function(x){
    dnbinom(as.numeric(names(obs)), mu = exp(x), size = theta_est)
  })
  nb_fit <- Reduce("+", density_list)
  vcd::rootogram(obs, nb_fit, main = paste0("Negative binomial fit for ", title_str), 
                 sub = paste0("Theta Estimate: ", round(theta_est, 2)))
}

plot_norm_counts <- function(exp_count_mat, gene_name, raw_var_estimate, title_str, 
                             scale_factor = scale_fac, xlim = c(0, 50)){
  par(mar = c(6,6,6,6))
  size_factors <- colSums(exp_count_mat) / scale_factor
  count_vals <- exp_count_mat[gene_name, ]
  norm_counts <- count_vals / size_factors
  lab_str <- paste0("mean: ", round(mean(norm_counts), 2), "\n",
                    "var: ", round(var(norm_counts), 2), "\n",
                    "est raw var: ", round(raw_var_estimate, 2))
  x <- hist(norm_counts, breaks = 50, main = paste0("Normalized Counts for ", title_str), xlim = xlim,
            cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
  text(xlim[2] * 0.5, y = max(x$counts) * 0.6, 
       labels = lab_str, pos = 4, cex = 1.5)
}

plot_pois_fit <- function(exp_count_mat, gene_name, title_str, scale_factor = scale_fac){
  size_factors <- colSums(exp_count_mat) / scale_factor
  count_vals <- exp_count_mat[gene_name, ]
  cs_mean <- mean(count_vals / size_factors)
  offs <- log(cs_mean * size_factors)
  obs <- table(count_vals)
  density_list <- lapply(offs, function(x){
    dpois(as.numeric(names(obs)), lambda = exp(x))
  })
  pois_fit <- Reduce("+", density_list)
  vcd::rootogram(obs, pois_fit, main = paste0("Poisson fit for ", title_str))
}

local_plot_ward_heatmap <- function(pp, annotation_col = NULL, annotation_row = NULL, 
                                    pal = NULL, cap_val = Inf){
  col_dist <- dist(t(pp), method = "euclidean")
  col_hclust <- hclust(col_dist, method = "ward.D2")
  row_dist <- dist(pp, method = "euclidean")
  row_hclust <- hclust(row_dist, method = "ward.D2")
  p <- pheatmap(pp, show_rownames = FALSE, annotation_col = annotation_col, 
                annotation_row = annotation_row,cluster_cols = col_hclust, 
                cluster_rows = row_hclust, 
                annotation_colors = pal)
}

################################################################################
# estimate SCV
################################################################################

## raw variance estimates
if(file.exists("R_Data/SC_SCV_Est.RData")){
  message("Loading data")
  obj_scv <- load("R_Data/SC_SCV_Est.RData")
} else {
  scv_est <- lapply(levels(scrna_sample_table$sample3), function(x){
    message(x)
    ids <- sc$sample3 == x
    sub <- sc[, ids]
    est_scv(sub[["RNA"]]@counts)
  })
  names(scv_est) <- levels(scrna_sample_table$sample3)
  save(scv_est, file = "R_Data/SC_SCV_Est.RData")
}

## select consistently highly expressed genes (min av. norm count > 1), excluding ribosomal proteins,
## and those not present in bulk RNA-seq and genes with raw variance estimate < 0
mean_mat <- sapply(scv_est, function(x) x$genes[, "cs_mean"])
raw_var_mat <- sapply(scv_est, function(x) x$genes[, "est_raw_var"])
min_vals <- apply(mean_mat, 1, min)
min_rv_vals <- apply(raw_var_mat, 1, min)
mean_vals <- rowMeans(mean_mat)
names(min_vals) <- rownames(sc)
names(min_rv_vals) <- rownames(sc)
names(mean_vals) <- rownames(sc)
keep <- min_vals > min_exp_val
table(keep)
gene_ids <- rownames(sc)[keep]
keep <- !grepl("^RP", gene_ids)
table(keep)
gene_ids <- gene_ids[keep]

## get inferred raw variance and mean values for these genes in long format
mean_exp <- sapply(scv_est, function(x) x$genes[gene_ids, "cs_mean"])
rownames(mean_exp) <- gene_ids
mean_vals <- rowMeans(mean_exp)
mean_df <- melt(mean_exp)
colnames(mean_df) <- c("gene", "sample", "mean_exp")

raw_var  <- sapply(scv_est, function(x) x$genes[gene_ids, "est_raw_var"])
rownames(raw_var) <- gene_ids
raw_var_df <- melt(raw_var)
colnames(raw_var_df) <- c("gene", "sample", "est_raw_var")

scv <- sapply(scv_est, function(x) x$genes[gene_ids, "est_scv"])
rownames(scv) <- gene_ids
scv_df <- melt(scv)
colnames(scv_df) <- c("gene", "sample", "est_scv")

xx_full <- cbind(mean_df, est_raw_var = raw_var_df$est_raw_var, est_scv = scv_df$est_scv)
type_ids <- match(xx_full$sample, scrna_sample_table$sample3)
xx_full$tnbc_type <- scrna_sample_table$tnbc.subtype[type_ids]
bin_vals <- 10^seq(floor(log10(min_exp_val)), 4, 1)
xx_full$mean_bin <- cut(mean_vals[as.character(xx_full$gene)], bin_vals)
xx_full$bas <- xx_full$gene %in% tnbc_gsl$bas$up
xx_full$lum <- xx_full$gene %in% tnbc_gsl$lum$up
xx_full$mes <- xx_full$gene %in% tnbc_gsl$mes$up
xx_full$de <- xx_full$gene %in% unlist(tnbc_gsl)
## filter out excluded (replicated) samples
keep <- !xx_full$sample %in% test_exc_samples
xx <- droplevels(xx_full[keep, ])

################################################################################
# analyse variance using summary stat method
################################################################################

## aggregate to sample level
agg_all <- aggregate(cbind(mean_exp, est_raw_var, est_scv) ~ sample,
                     data = xx, FUN = mean, drop = FALSE)
agg_non_de <- aggregate(cbind(mean_exp, est_raw_var, est_scv) ~ sample, 
                        data = subset(xx, !de), FUN = mean, drop = FALSE)
agg_mes <- aggregate(cbind(mean_exp, est_raw_var, est_scv) ~ sample,
                     data = subset(xx, mes), FUN = mean, drop = FALSE)
agg_lum <- aggregate(cbind(mean_exp, est_raw_var, est_scv) ~ sample,
                     data = subset(xx, lum), FUN = mean, drop = FALSE)

type_ids <- match(agg_non_de$sample, scrna_sample_table$sample3)
agg_all$tnbc_type <- factor(scrna_sample_table$tnbc.subtype[type_ids], levels = c("basal", "luminal", "mesenchymal"))
agg_non_de$tnbc_type <- factor(scrna_sample_table$tnbc.subtype[type_ids], levels = c("basal", "luminal", "mesenchymal"))
agg_mes$tnbc_type <- factor(scrna_sample_table$tnbc.subtype[type_ids], levels = c("basal", "luminal", "mesenchymal"))
agg_lum$tnbc_type <- factor(scrna_sample_table$tnbc.subtype[type_ids], levels = c("basal", "luminal", "mesenchymal"))

res_all <- list()
pdf("analysis/scrna_seq/diversity/all_boxplot_mean_expression.pdf")
res_all[["mean"]] <- local_boxplot(agg_all, "mean_exp", "All Genes Mean Expression", "Av. Est. Mean Expression", log_y = TRUE)
dev.off()
pdf("analysis/scrna_seq/diversity/all_boxplot_raw_var.pdf")
res_all[["rv"]] <- local_boxplot(agg_all, "est_raw_var", "All Genes Raw Var", "Mean Est. Raw Var")
dev.off()
pdf("analysis/scrna_seq/diversity/all_boxplot_scv.pdf")
res_all[["scv"]] <- local_boxplot(agg_all, "est_scv", "All Genes Raw SCV", "Mean Est. Raw SCV")
dev.off()

res_non_de <- list()
pdf("analysis/scrna_seq/diversity/non_de_boxplot_mean_expression.pdf")
res_non_de[["mean"]] <- local_boxplot(agg_non_de, "mean_exp", "Non-DE Genes Mean Expression", "Av. Est. Mean Expression", log_y = TRUE)
dev.off()
pdf("analysis/scrna_seq/diversity/non_de_boxplot_raw_var.pdf")
res_non_de[["rv"]] <- local_boxplot(agg_non_de, "est_raw_var", "Non-DE Genes Raw Var", "Mean Est. Raw Var")
dev.off()
pdf("analysis/scrna_seq/diversity/non_de_boxplot_scv.pdf")
res_non_de[["scv"]] <- local_boxplot(agg_non_de, "est_scv", "Non-DE Genes Raw SCV", "Mean Est. Raw SCV")
dev.off()

res_lum <- list()
pdf("analysis/scrna_seq/diversity/lum_up_boxplot_mean_expression.pdf")
res_lum[["mean"]] <- local_boxplot(agg_lum, "mean_exp", "Lum Genes Mean Expression", "Av. Est. Mean Expression", log_y = TRUE)
dev.off()
pdf("analysis/scrna_seq/diversity/lum_up_boxplot_raw_var.pdf")
res_lum[["rv"]] <- local_boxplot(agg_lum, "est_raw_var", "Lum Genes Raw Var", "Mean Est. Raw Var")
dev.off()
pdf("analysis/scrna_seq/diversity/lum_up_boxplot_scv.pdf")
res_lum[["scv"]] <- local_boxplot(agg_lum, "est_scv", "Lum Genes Raw SCV", "Mean Est. Raw SCV")
dev.off()

res_mes <- list()
pdf("analysis/scrna_seq/diversity/mes_up_boxplot_mean_expression.pdf")
res_mes[["mean"]] <- local_boxplot(agg_mes, "mean_exp", "Mes Genes Mean Expression", "Av. Est. Mean Expression", log_y = TRUE)
dev.off()
pdf("analysis/scrna_seq/diversity/mes_up_boxplot_raw_var.pdf")
res_mes[["rv"]] <- local_boxplot(agg_mes, "est_raw_var", "Mes Genes Raw Var", "Mean Est. Raw Var")
dev.off()
pdf("analysis/scrna_seq/diversity/mes_up_boxplot_scv.pdf")
res_mes[["scv"]] <- local_boxplot(agg_mes, "est_scv", "Mes Genes Raw SCV", "Mean Est. Raw SCV")
dev.off()

bin_agg_all <- aggregate(cbind(mean_exp, est_raw_var, est_scv) ~ sample + mean_bin,
                         data = xx, FUN = mean, drop = FALSE)
bin_agg_non_de <- aggregate(cbind(mean_exp, est_raw_var, est_scv) ~ sample + mean_bin, 
                                  data = subset(xx, !de), FUN = mean, drop = FALSE)
bin_agg_mes <- aggregate(cbind(mean_exp, est_raw_var, est_scv) ~ sample + mean_bin,
                               data = subset(xx, mes), FUN = mean, drop = FALSE)
bin_agg_lum <- aggregate(cbind(mean_exp, est_raw_var, est_scv) ~ sample + mean_bin,
                               data = subset(xx, lum), FUN = mean, drop = FALSE)

type_ids <- match(bin_agg_non_de$sample, scrna_sample_table$sample3)
bin_agg_all$tnbc_type <- factor(scrna_sample_table$tnbc.subtype[type_ids], levels = c("basal", "luminal", "mesenchymal"))
bin_agg_non_de$tnbc_type <- factor(scrna_sample_table$tnbc.subtype[type_ids], levels = c("basal", "luminal", "mesenchymal"))
bin_agg_mes$tnbc_type <- factor(scrna_sample_table$tnbc.subtype[type_ids], levels = c("basal", "luminal", "mesenchymal"))
bin_agg_lum$tnbc_type <- factor(scrna_sample_table$tnbc.subtype[type_ids], levels = c("basal", "luminal", "mesenchymal"))

res_bin_all <- list()
pdf("analysis/scrna_seq/diversity/all_bin_boxplot_mean_expression.pdf")
res_bin_all[["mean"]] <- local_facet_boxplot(bin_agg_all, "mean_exp", "All Genes Mean Expression", "Av. Est. Mean Expression", log_y = TRUE)
dev.off()
pdf("analysis/scrna_seq/diversity/all_bin_boxplot_raw_var.pdf")
res_bin_all[["rv"]] <- local_facet_boxplot(bin_agg_all, "est_raw_var", "All Genes Raw Var", "Mean Est. Raw Var")
dev.off()
pdf("analysis/scrna_seq/diversity/all_bin_boxplot_scv.pdf")
res_bin_all[["scv"]] <- local_facet_boxplot(bin_agg_all, "est_scv", "All Genes Raw SCV", "Mean Est. Raw SCV")
dev.off()

res_bin_non_de <- list()
pdf("analysis/scrna_seq/diversity/non_de_bin_boxplot_mean_expression.pdf")
res_bin_non_de[["mean"]] <- local_facet_boxplot(bin_agg_non_de, "mean_exp", "Non-DE Genes Mean Expression", "Av. Est. Mean Expression", log_y = TRUE)
dev.off()
pdf("analysis/scrna_seq/diversity/non_de_bin_boxplot_raw_var.pdf")
res_bin_non_de[["rv"]] <- local_facet_boxplot(bin_agg_non_de, "est_raw_var", "Non-DE Genes Raw Var", "Mean Est. Raw Var")
dev.off()
pdf("analysis/scrna_seq/diversity/non_de_bin_boxplot_scv.pdf")
res_bin_non_de[["scv"]] <- local_facet_boxplot(bin_agg_non_de, "est_scv", "Non-DE Genes Raw SCV", "Mean Est. Raw SCV")
dev.off()

res_bin_lum <- list()
pdf("analysis/scrna_seq/diversity/lum_up_bin_boxplot_mean_expression.pdf")
res_bin_lum[["mean"]] <- local_facet_boxplot(bin_agg_lum, "mean_exp", "Lum Genes Mean Expression", "Av. Est. Mean Expression", log_y = TRUE)
dev.off()
pdf("analysis/scrna_seq/diversity/lum_up_bin_boxplot_raw_var.pdf")
res_bin_lum[["rv"]] <- local_facet_boxplot(bin_agg_lum, "est_raw_var", "Lum Genes Raw Var", "Mean Est. Raw Var")
dev.off()
pdf("analysis/scrna_seq/diversity/lum_up_bin_boxplot_scv.pdf")
res_bin_lum[["scv"]] <- local_facet_boxplot(bin_agg_lum, "est_scv", "Lum Genes Raw SCV", "Mean Est. Raw SCV")
dev.off()

res_bin_mes <- list()
pdf("analysis/scrna_seq/diversity/mes_up_bin_boxplot_mean_expression.pdf")
res_bin_mes[["mean"]] <- local_facet_boxplot(bin_agg_mes, "mean_exp", "Mes Genes Mean Expression", "Av. Est. Mean Expression", log_y = TRUE)
dev.off()
pdf("analysis/scrna_seq/diversity/mes_up_bin_boxplot_raw_var.pdf")
res_bin_mes[["rv"]] <- local_facet_boxplot(bin_agg_mes,"est_raw_var", "Mes Genes Raw Var", "Mean Est. Raw Var")
dev.off()
pdf("analysis/scrna_seq/diversity/mes_up_bin_boxplot_scv.pdf")
res_bin_mes[["scv"]] <- local_facet_boxplot(bin_agg_mes, "est_scv", "Mes Genes Raw SCV", "Mean Est. Raw SCV")
dev.off()

save(res_all,res_non_de, res_lum, res_mes, res_bin_non_de, res_bin_lum, res_bin_mes, 
     file = "R_Data/SC_Diversity_Comparisons.RData")

################################################################################
# compare estimates for replicates
################################################################################

col_ann <- scrna_sample_table[, "tnbc.subtype", drop = FALSE]
rownames(col_ann) <- scrna_sample_table$sample3
row_ann <- data.frame(row.names = rownames(scv))
lab_genes <- c("UQCRH", "HIST1H4C")
row_ann$gene <- factor(sapply(rownames(row_ann), function(x){
  if(x %in% lab_genes){
    x
  } else {
    "Other"
  }
}), levels = c(lab_genes, "Other"))
gene_cols <- c("blue", "red", "grey")
names(gene_cols) <- c(lab_genes, "Other")
pal <- list(tnbc.subtype = line_colours$type, gene = gene_cols)

pdf("analysis/scrna_seq/diversity/heatmap_ward_mean.pdf")
p <- local_plot_ward_heatmap(log(mean_exp), annotation_col = col_ann, pal = pal)
print(p)
dev.off()

pdf("analysis/scrna_seq/diversity/heatmap_ward_scv.pdf")
p <- local_plot_ward_heatmap(scv, annotation_col = col_ann, pal = pal)
print(p)
dev.off()

pdf("analysis/scrna_seq/diversity/heatmap_asinh_scv.pdf")
p <- local_plot_ward_heatmap(asinh(scv * 5), annotation_col = col_ann, 
                             annotation_row = row_ann, pal = pal)
print(p)
dev.off()

pdf("analysis/scrna_seq/diversity/heatmap_asinh_raw_var.pdf")
p <- local_plot_ward_heatmap(asinh(raw_var * 5), annotation_col = col_ann, 
                             annotation_row = row_ann, pal = pal)
print(p)
dev.off()
################################################################################
# analysis based on theta estimates from negative binomial regression
################################################################################

## Analyze theta estimates based on negative binomial regression
if(file.exists("R_Data/SC_NB_Reg.RData")){
  message("Loading data")
  load("R_Data/SC_NB_Reg.RData")
} else {
  nb_est <- lapply(levels(scrna_sample_table$sample3), function(x){
    message(x)
    ids <- sc$sample3 == x
    sub <- sc[, ids]
    est_nb_reg(sub[["RNA"]]@counts, gene_ids)
  })
  names(nb_est) <- levels(scrna_sample_table$sample3)
  save(nb_est, file = "R_Data/SC_NB_Reg.RData")
}

nb_theta_ests <- sapply(nb_est, function(x) x$genes$est_theta)
rownames(nb_theta_ests) <- gene_ids
## remove genes where with regression warnings
keep <- !colnames(nb_theta_ests) %in% test_exc_samples
nb_theta_ests_test <- nb_theta_ests[, keep]
keep <- apply(nb_theta_ests_test, 1, function(x) length(which(is.na(x))) == 0)
nb_theta_ests_narm <- nb_theta_ests_test[keep, ]

nb_theta_narm_df <- melt(nb_theta_ests_narm)
colnames(nb_theta_narm_df) <- c("gene", "sample", "theta_est")
type_ids <- match(nb_theta_narm_df$sample, scrna_sample_table$sample3)
nb_theta_narm_df$tnbc_type <- scrna_sample_table$tnbc.subtype[type_ids]
nb_theta_narm_df$de <- nb_theta_narm_df$gene %in% unlist(tnbc_gsl)
nb_theta_narm_df$bas <- nb_theta_narm_df$gene %in% tnbc_gsl$bas$up
nb_theta_narm_df$lum <- nb_theta_narm_df$gene %in% tnbc_gsl$lum$up
nb_theta_narm_df$mes <- nb_theta_narm_df$gene %in% tnbc_gsl$mes$up

#agg_non_de <- aggregate(cbind(mean_exp, est_raw_var, est_scv) ~ sample, 
#                        data = subset(xx, !de), FUN = mean, drop = FALSE)

#nb_theta_narm_agg <- aggregate(theta_est ~ sample, data = subset(nb_theta_narm_df, !de), FUN = mean)
nb_theta_narm_agg <- aggregate(theta_est ~ sample, data = nb_theta_narm_df, FUN = mean)
type_ids <- match(nb_theta_narm_agg$sample, scrna_sample_table$sample3)
nb_theta_narm_agg$tnbc_type <- factor(scrna_sample_table$tnbc.subtype[type_ids], 
                                      levels = c("basal", "luminal", "mesenchymal"))

pdf("analysis/scrna_seq/diversity/all_boxplot_theta.pdf")
res_theta <- local_boxplot(nb_theta_narm_agg, "theta_est", "Non-DE genes", "Est. Theta", 
                           cap_val = 1e6)
dev.off()

## compare to mom method
nb_theta_df <- melt(nb_theta_ests)
colnames(nb_theta_df) <- c("gene", "sample", "est_theta")

pp <- data.frame(mom_theta_est = 1 / scv_df$est_scv, nb_theta_est = nb_theta_df$est_theta)
pp$mom_theta_est <- pmin(pp$mom_theta_est, 100)
pp$nb_theta_est <- pmin(pp$nb_theta_est, 100)
pdf("analysis/scrna_seq/diversity/theta_est_comparison.pdf")
plot(pp$mom_theta_est, pp$nb_theta_est, xlab = "Method of Moments Inverse Raw Squared CV", ylab = "NB Regression Theta",
     main = "Variability Estaimtes")
abline(a = 0, b = 1)
dev.off()

nb_theta_df$cs_mean <- mean_df$mean_exp
sort_ids <- intersect(order(nb_theta_df$est_theta, decreasing = TRUE), which(nb_theta_df$cs_mean > 10))
sort_ids_non_na <- intersect(sort_ids, which(!is.na(nb_theta_df$est_theta)))
(hi_ex <- nb_theta_df[head(sort_ids_non_na, 1), ])
(lo_ex <- nb_theta_df[tail(sort_ids_non_na, 1), ])
(hi_raw_var_est <- scv_est[[as.character(hi_ex$sample)]]$genes[as.character(hi_ex$gene), "est_raw_var"])
(lo_raw_var_est <- scv_est[[as.character(lo_ex$sample)]]$genes[as.character(lo_ex$gene), "est_raw_var"])

sub <- sc[, sc$sample3 == hi_ex$sample]
pdf("analysis/scrna_seq/diversity/hi_norm_counts_example.pdf")
plot_norm_counts(sub[["RNA"]]@counts, hi_ex$gene, hi_raw_var_est, paste0(hi_ex$gene, " in ", hi_ex$sample))
dev.off()
pdf("analysis/scrna_seq/diversity/hi_nbin_example.pdf")
plot_nb_fit(sub[["RNA"]]@counts, hi_ex$gene, paste0(hi_ex$gene, " in ", hi_ex$sample))
dev.off()
pdf("analysis/scrna_seq/diversity/hi_pois_example.pdf")
plot_pois_fit(sub[["RNA"]]@counts, hi_ex$gene, paste0(hi_ex$gene, " in ", hi_ex$sample))
dev.off()

sub <- sc[, sc$sample3 == lo_ex$sample]
pdf("analysis/scrna_seq/diversity/lo_norm_counts_example.pdf")
plot_norm_counts(sub[["RNA"]]@counts, lo_ex$gene, lo_raw_var_est, paste0(lo_ex$gene, " in ", lo_ex$sample))
dev.off()
pdf("analysis/scrna_seq/diversity/lo_nbin_example.pdf")
plot_nb_fit(sub[["RNA"]]@counts, lo_ex$gene, paste0(lo_ex$gene, " in ", lo_ex$sample))
dev.off()
pdf("analysis/scrna_seq/diversity/lo_pois_example.pdf")
plot_pois_fit(sub[["RNA"]]@counts, lo_ex$gene, paste0(lo_ex$gene, " in ", lo_ex$sample))
dev.off()

################################################################################
# analyse variance using mixed models
################################################################################

min_raw_var <- aggregate(est_raw_var ~ gene, xx, min)
keep_genes <- min_raw_var$gene[which(min_raw_var$est_raw_var > 0)]
xx_m <- droplevels(subset(xx, gene %in% keep_genes))

## individual gene regressions
lm_coeffs <- as.data.frame(t(sapply(levels(xx_m$gene), function(x){
  xx_m_sub <- subset(xx_m, gene == x)
  #lm_res <- lm(log(est_raw_var + 1) ~ log(mean_exp + 1), xx_m_sub)
  lm_res <- lm(log(est_raw_var) ~ log(mean_exp), xx_m_sub)
  lm_res$coefficients
})))
colnames(lm_coeffs) <- c("intercept", "coefficient")
lm_coeffs$gene <- rownames(lm_coeffs)

## plot individual fits for panel of random genes
set.seed(1234)
rand_genes <- as.character(sample(unique(xx_m$gene), 10))

gg <- subset(xx_m, gene %in% rand_genes)

pdf("analysis/scrna_seq/diversity/reg_sample_genes.pdf")
ggplot(gg, aes(x = log(mean_exp), y = log(est_raw_var), color = gene)) +
  geom_point() +
  my_theme() +
  ggtitle("Random Genes")
dev.off()

pdf("analysis/scrna_seq/diversity/reg_sample_genes_enlarged.pdf", 12, 9)
ggplot(gg, aes(x = log(mean_exp), y = log(est_raw_var), color = gene)) +
  geom_point() +
  my_theme() +
  ggtitle("Random Genes")
dev.off()

pdf("analysis/scrna_seq/diversity/reg_sample_all_genes.pdf")
ggplot(xx_m, aes(x = log(mean_exp), y = log(est_raw_var), color = tnbc_type)) +
  geom_point(size = 0.1) +
  my_theme() +
  ggtitle("All Filtered Genes")
dev.off()


pdf("analysis/scrna_seq/diversity/reg_sample_fits.pdf")
par(mfrow = c(4, 3))
lapply(rand_genes, function(x){
  xx_m_sub <- subset(xx_m, gene == x)
  lm_sum <- summary(lm(log(est_raw_var) ~ log(mean_exp), xx_m_sub))
  plot(log(est_raw_var) ~ log(mean_exp), xx_m_sub, 
       col = line_colours$type[xx_m_sub$tnbc_type],
       #xlim = c(0, 6), ylim = c(-8, 10),
       main = x, xlab = "log(Mean Expression)", ylab = "log(Est. Raw Var)")
  abline(a = lm_sum$coefficients[1,1], b = lm_sum$coefficients[2,1])
})
dev.off()

pdf("analysis/scrna_seq/diversity/reg_individual_intercepts.pdf")
hist(lm_coeffs$intercept, xlab = "Intercept", main = "Regression Intercepts")
abline(v = mean(lm_coeffs$intercept))
dev.off()

pdf("analysis/scrna_seq/diversity/reg_individual_coeffs.pdf")
hist(lm_coeffs$coefficient, xlab = "log(Mean Expression) Coefficient", 
     main = "Regression Coefficients")
abline(v = mean(lm_coeffs$coefficient))
dev.off()

## Ordinary mixed model
#lmer_full_res <- lmer(log(est_raw_var) ~ log(mean_exp * 5) + tnbc_type + gene + (1|sample), 
#                      xx_m)
lmer_res <- lmer(log(est_raw_var) ~ log(mean_exp) + tnbc_type + (1|gene) + (1|sample), 
                 xx_m)
(lmer_anova <- Anova(lmer_res))
lmer_sum <- summary(lmer_res)
lmer_resid <- residuals(lmer_res)
(lmer_ci <- confint(lmer_res, parm = c("tnbc_typeluminal", "tnbc_typemesenchymal"), level = 0.95))

## forest plot
forest_data <- cbind(lmer_sum$coefficients[rownames(lmer_ci), 1], lmer_ci)
tabletext <- matrix(c("Luminal", "Mesenchymal"))

pdf("analysis/scrna_seq/diversity/reg_lmer_forestplot.pdf")
forestplot(tabletext, forest_data, lwd.zero = 2, lwd.ci = 2,
           txt_gp = fpTxtGp(ticks = gpar(cex = 2), 
                            label = gpar(cex = 2)))
dev.off()

#forestplot(tabletext, exp_forest_data, xlog = TRUE, 
#           txt_gp = fpTxtGp(ticks = gpar(cex = cex_tick), 
#                            label = gpar(cex = cex_lab)))

pdf("analysis/scrna_seq/diversity/reg_lmer_diagnostics.pdf")
plot(lmer_res)
hist(lmer_resid, xlab = "Residual", main = "Residuals")
dev.off()

lmer_nde_res <- lmer(log(est_raw_var) ~ log(mean_exp) + tnbc_type + (1|gene) + (1|sample), 
                     subset(xx_m, !de))
(lmer_nde_anova <- Anova(lmer_nde_res))
lmer_nde_sum <- summary(lmer_nde_res)
(lmer_nde_ci <- confint(lmer_nde_res, parm = c("tnbc_typeluminal", "tnbc_typemesenchymal"), level = 0.95))

lmer_de_res <- lmer(log(est_raw_var) ~ log(mean_exp) + tnbc_type + (1|gene) + (1|sample), 
                    subset(xx_m, de))
(lmer_de_anova <- Anova(lmer_de_res))
lmer_de_sum <- summary(lmer_de_res)
(lmer_de_ci <- confint(lmer_de_res, parm = c("tnbc_typeluminal", "tnbc_typemesenchymal"), level = 0.95))

lmer_lum_res <- lmer(log(est_raw_var) ~ log(mean_exp) + tnbc_type + (1|gene) + (1|sample), 
                     subset(xx_m, lum))
(lmer_lum_anova <- Anova(lmer_lum_res))
lmer_lum_sum <- summary(lmer_lum_res)
(lmer_lum_ci <- confint(lmer_lum_res, parm = c("tnbc_typeluminal", "tnbc_typemesenchymal"), level = 0.95))

lmer_mes_res <- lmer(log(est_raw_var) ~ log(mean_exp) + tnbc_type + (1|gene) + (1|sample), 
                     subset(xx_m, mes))
(lmer_mes_anova <- Anova(lmer_mes_res))
lmer_mes_sum <- summary(lmer_mes_res)
(lmer_mes_ci <- confint(lmer_mes_res, parm = c("tnbc_typeluminal", "tnbc_typemesenchymal"), level = 0.95))

save(lmer_res, lmer_ci, lmer_de_res, lmer_de_ci, file = "R_Data/SC_LMER.RData")

## random coefficients model
lmer_rcoef_res <- lmer(log(est_raw_var) ~ log(mean_exp) + tnbc_type + (log(mean_exp) +1|gene) + (log(mean_exp)|sample), 
                       xx_m)
(lmer_rcoef_anova <- Anova(lmer_rcoef_res))
lmer_rcoef_sum <- summary(lmer_rcoef_res)
lmer_rcoef_resid <- residuals(lmer_rcoef_res)
(lmer_rcoef_ci <- confint(lmer_rcoef_res, parm = c("tnbc_typeluminal", "tnbc_typemesenchymal"), level = 0.95))

## GLMM
glmer_res <- glmer(est_raw_var ~ log(mean_exp) + tnbc_type + (1|gene) + (1|sample),
                   xx_m, family = Gamma(link = "log"))
glmer_resid <- residuals(glmer_res, type = "pearson", scaled = TRUE)
(glmer_sum <- summary(glmer_res))



