source("../utilities/utils.R")
library(survminer)
library(survival)
library(ggplot2)
library(pheatmap)
library(chisq.posthoc.test)
library(FSA)

load("../R_Data/Heatmap_Metadata.RData")

#----------------------------------------------------------------------------------------
#' Check procuerment status for BC sub-type for TCGA
#' data
tcga_log_test <- function(val){
  test_stat <- xmlAttrs(val)["procurement_status"]
  if(test_stat == "Completed"){
    res <- xmlValue(val)
  } else {
    res <- test_stat
  }
  return(res)
}

#----------------------------------------------------------------------------------------
#' Classify samples as tnbc, non-tnbc or indeterminate
classify_tnbc <- function(er, pr, her2){
  if ((er == "Negative") & (pr == "Negative") & (her2 == "Negative")){
    res <- "tnbc"
  } else if ((er=="Positive") | (pr == "Positive") | (her2 == "Positive")) {
    res <- "non-tnbc"
  } else {
    res <- "indeterminate"
  }
  return(res)
}

#----------------------------------------------------------------------------------------
#' Plot violin
plot_vioplot <- function(type, y, ylab, ann_col = line_colours$type2, text_size = 50, lab_size = 10){
  stopifnot(is.factor(type))
  
  data <- data.frame(type = type, y = as.numeric(y))
  data_summary <- data.frame(type = levels(data$type))
  data_summary$x <- 1:nrow(data_summary)
  data_summary$n <- paste0("N=", sapply(data_summary$type, function(x) length(which(data$type == x))))
  data_summary$max.y <- sapply(data_summary$type, function(x) max(data$y[which(data$type == x)]))
  
  ggplot(data, aes(x = type, y = y, fill = type)) +
    geom_violin() +
    scale_fill_manual(values = ann_col) +
    ylab(ylab) +
    theme(text = element_text(size = text_size)) +
    geom_text(size = lab_size,
              data = data_summary, 
              show.legend = FALSE, 
              aes(x = x, y = max.y * 1.2, label = n))
}

#----------------------------------------------------------------------------------------
#' Plot boxplot
plot_boxplot <- function(type, y, ylab, ann_col = line_colours$type2, text_size = 50, lab_size = 10){
  stopifnot(is.factor(type))
  
  data <- data.frame(type = type, y = as.numeric(y))
  data_summary <- data.frame(type = levels(data$type))
  data_summary$x <- 1:nrow(data_summary)
  data_summary$n <- paste0("N=", sapply(data_summary$type, function(x) length(which(data$type == x))))
  data_summary$max.y <- sapply(data_summary$type, function(x) max(data$y[which(data$type == x)]))
  
  ggplot(data, aes(x = type, y = y, fill = type)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = ann_col) +
    geom_jitter(color = "black", size = 3, alpha = 0.9) +
    theme(text = element_text(size = text_size)) +
    ylab(ylab) +
    geom_text(size = lab_size,
              data = data_summary,
              show.legend = FALSE,
              aes(x = x, y = max.y * 1.2, label = n))
}

#----------------------------------------------------------------------------------------
#' plot scatter plot for y vs. trait value
plot_scatter <- function(trait, y, xlab, ylab, text_size = 60){
  data <- data.frame(x = trait, y = y)
  ggplot(data, aes(x = trait, y = y)) +
    geom_point() +
    xlab(xlab) +
    ylab(ylab) +
    theme(text = element_text(size = text_size))
}

#----------------------------------------------------------------------------------------
#' km_plot
#' 
#' Generate Kaplan-Meier plot
km_plot <- function(status_id, km_y_lab, time_id, var_id, f_surv_data, km_colors = NULL, 
                    km_x_lab = "Time (Months)", fs = 3, p_val = TRUE, p_val_method = TRUE, 
                    text_size = 10){
  # asserts
  stopifnot(is.null(km_colors) | is.character(km_colors))
  validate_survival_data(status_id, time_id, var_id, f_surv_data)
  # plot
  if(!is.null(km_colors)) {
    stopifnot(all(levels(f_surv_data[[var_id]]) %in% names(km_colors)))
    km_colors <- as.character(km_colors[levels(f_surv_data[[var_id]])])
  }
  fml <- as.formula(paste0("Surv(", time_id, ", ", status_id, ") ~ ", var_id))
  fit1 <- surv_fit(fml, f_surv_data)
  #tab_height <- length(levels(f_surv_data[[var_id]])) * 0.1 + 0.1
  tab_height = 0.3
  strata_names <- sapply(names(fit1$strata), function(x) strsplit(x, "=")[[1]][2])
  stopifnot(all(strata_names == levels(f_surv_data[[var_id]])))
  p <- ggsurvplot(fit = fit1, data = f_surv_data, pval = p_val, pval.method = p_val_method, risk.table = TRUE, 
                  ggtheme = my_theme(text_size = text_size), xlab = km_x_lab, ylab = km_y_lab, pval.size = fs, 
                  fontsize = fs, palette = km_colors, legend.labs = levels(f_surv_data[[var_id]]),
                  tables.height = tab_height)
  return(p)
}

#----------------------------------------------------------------------------------------
#' cox_regression
#' 
#' Run univariate or multivariate cox regression, with or without covariates
cox_regression <- function(status_id, time_id, var_ids, f_surv_data, cov_ids = c()){
  validate_survival_data(status_id, time_id, var_ids, f_surv_data, cov_ids)
  fml <- paste0("Surv(", time_id, ", ", status_id, ") ~ ", paste(c(var_ids, cov_ids), collapse = " + "))
  cox_out <- coxph_summary(coxph(as.formula(fml), f_surv_data))
  return(cox_out)
}

#----------------------------------------------------------------------------------------
#' write_kruskal_test
#' 
#' Output kruskal-Wallis test results to file
write_kruskal_test <- function(type, y, path, comp_id, alpha = 0.05){
  df <- data.frame(x = as.numeric(y), g = type)
  kst <- kruskal.test(x ~ g, data = df)$p.value
  write.table(kst, file = paste0(path, "/", comp_id, "_kruskal_test_p_value.csv"), sep = ",")
  agg_res <- aggregate(x ~ g, data = df, FUN = median)
  write.table(agg_res, file = paste0(path, "/", comp_id, "_group_medians.csv"), sep = ",", col.names = NA)
  if(kst < alpha){
    d_res <- dunnTest(x ~ g, data = df, method = "holm")$res
    write.table(d_res, file = paste0(path, "/", comp_id, "_dunn_test.csv"), sep = ",", col.names = NA)  
  }
}

#----------------------------------------------------------------------------------------
#' write_wilcox_test
#' 
#' Output Wilcox test results to file
write_wilcox_test <- function(type, y, path, comp_id){
  df <- data.frame(x = as.numeric(y), g = type)
  wc_res <- aggregate(x ~ g, data = df, FUN = median)
  write.table(wc_res, file = paste0(path, "/", comp_id, "_group_medians.csv"), sep = ",", col.names = NA)
  wct <- wilcox.test(x ~ g, data = df)$p.value
  write.table(wct, file = paste0(path, "/", comp_id, "_wilcox_test_p_value.csv"))
}

#----------------------------------------------------------------------------------------
#' assess_gsl_regulation
#' 
#' Generate vioplots of absolute correlations and raw pearson correlations with
#' gene expression for candidate up and down targets in each candidate target
#' set. Test for difference between values for up and down targets using
#' Wilcoxon Rank sum test in each case
assess_gsl_regulation <- function(exp_data, regulator, gsl, plot_path, list_path){
  # check variance of regulator is non-zero for correlation analysis
  reg_exp <- as.numeric(exp_data[regulator, ])
  stopifnot(var(reg_exp) > 0)
  # only consider correlations with genes that have non-zero variance
  keep <- apply(exp_data, 1, var) > 0
  n_up_dn_sig_keep <- sapply(gsl[[1]], function(x) length(which(rownames(exp_data)[keep] %in% x)))
  stopifnot(all(n_up_dn_sig_keep > 0))
  exp_data <- exp_data[keep, ]
  
  reg_cors <- data.frame(row.names = rownames(exp_data), cor = as.numeric(cor(t(exp_data), reg_exp)))
  reg_cors$abs.cor <- abs(reg_cors$cor)
  
  reg_cors$signature <- rownames(reg_cors) %in% unlist(gsl[[1]])
  
  pdf(paste0(plot_path, names(gsl), "_abs_cor.pdf"), width = 11, height = 9)
  p <- plot_vioplot(type = factor(reg_cors$signature, levels = c(FALSE, TRUE)), 
                    y = reg_cors$abs.cor, ann_col = c("blue", "red"), ylab = "Abs Correlation")
  print(p)
  dev.off()
  
  sig_reg_cors <- reg_cors[which(reg_cors$signature),]
  sig_reg_cors$dir <- sapply(rownames(sig_reg_cors), function(x) {
    if(x %in% gsl[[1]]$up) {
      "up"
    } else {
      "dn"
    }
  })
  
  pdf(paste0(plot_path, names(gsl), "_cor.pdf"), width = 11, height = 9)
  p <- plot_vioplot(type = factor(sig_reg_cors$dir, levels = c("dn", "up")), 
                    y = sig_reg_cors$cor, ann_col = c("blue", "red"), ylab = "Correlation")
  print(p)
  dev.off()
  
  write_wilcox_test(type = reg_cors$signature, y = reg_cors$abs.cor, path = list_path, comp_id = paste0(names(gsl), "_abs_cor"))
  write_wilcox_test(type = sig_reg_cors$dir, y = sig_reg_cors$cor, path = list_path, comp_id = paste0(names(gsl), "_cor"))
}

#----------------------------------------------------------------------------------------
#' assess_gml_regulation
#' 
#' Calculate absolute correlations for each gene in each candidate target module, and plot
#' with vioplot vs. absolute correlations for other genes. Test for differences using
#' Wilcoxon Rank Sum test
assess_gml_regulation <- function(exp_data, regulator, gml, plot_path, list_path){
  reg_exp <- as.numeric(exp_data[regulator, ])
  stopifnot(var(reg_exp) > 0)
  keep <- apply(exp_data, 1, var) > 0
  n_module_keep <- length(which(rownames(exp_data)[keep] %in% gml[[1]]))
  stopifnot(n_module_keep > 0)
  exp_data <- exp_data[keep, ]
  
  reg_cors <- data.frame(row.names = rownames(exp_data), abs.cor = abs(as.numeric(cor(t(exp_data), reg_exp))))
  
  reg_cors$signature <- rownames(reg_cors) %in% unlist(gml)
  
  pdf(paste0(plot_path, names(gml), "_abs_cor.pdf"), width = 11, height = 9)
  p <- plot_vioplot(type = factor(reg_cors$signature, levels = c(FALSE, TRUE)),
                    y = reg_cors$abs.cor, ann_col = c("blue", "red"), ylab = "Abs Correlation")
  print(p)
  dev.off()
  
  write_wilcox_test(type = reg_cors$signature, y = reg_cors$abs.cor, path = list_path, comp_id = paste0(names(gml), "_abs_cor"))
}

#----------------------------------------------------------------------------------------
#' analyse_signature_by_type
#' 
#' Plots vioplot of signature score against type for each signature. Tests for 
#' associations using Kruskall-Wallis tests. Associations are adjusted by Holm's Post-hoc 
#' analysis using Dunn Test is carried out for the significant associations before filtering
analyse_signature_by_type <- function(type, sig_scores, path, category_id, adjust_p_vals = TRUE, 
                                      alpha = 0.05, ann_col = line_colours$type2, text_size = 50,
                                      lab_size = 10, sig_meta = NULL, n_sum_col = 2){
  plot_path <- paste0(path, "plots/")
  list_path <- paste0(path, "lists/")
  summary_path <- paste0(path, "summary/")
  if(!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
  if(!dir.exists(list_path)) dir.create(list_path, recursive = TRUE)
  if(!dir.exists(summary_path)) dir.create(summary_path, recursive = TRUE)
  
  for(i in 1:ncol(sig_scores)){
    sig_id <- colnames(sig_scores)[i]
    message(sig_id)
    
    pdf(paste0(plot_path, sig_id, "_by_", category_id, ".pdf"), width = 11, height = 9)
    p <- plot_vioplot(type = type, y = sig_scores[,i], ylab = paste0(sig_id, " Score"), 
                      ann_col = ann_col, text_size = text_size, lab_size = lab_size)
    print(p)
    dev.off()
  }
  
  stat_vals <- apply(sig_scores, 2, function(x){
    kruskal.test(x, type)$statistic
  })
  p_vals <- apply(sig_scores, 2, function(x){
    kruskal.test(x, type)$p.value
  })
  sig_tests <- data.frame(stat = stat_vals, p.val = p_vals)
  if(adjust_p_vals){
    sig_tests$p.adj <- p.adjust(p_vals, method = "holm")
  } 
  
  sig_tests <- sig_tests[order(sig_tests$stat, decreasing = TRUE), ]
  
  write.table(sig_tests, file = paste0(summary_path, category_id, "_kruskal_test_p_values.csv"), sep = ",", col.names = NA)
  
  signif_tests <- sig_tests[which(sig_tests$p.val < alpha), ]
  
  if(nrow(signif_tests) > 0){
    for(i in 1:nrow(signif_tests)){
      sig_id <- rownames(signif_tests)[i]
      sig_vals <- sig_scores[, match(sig_id, colnames(sig_scores))]
      write_kruskal_test(type = type, y = sig_vals, path = list_path, comp_id = paste0(sig_id, "_by_", category_id))
    }  
  }
  
  # create summary heatmap
  aa <- aggregate(sig_scores, by = list(type), median)
  at <- t(aa[,-1, drop = FALSE])
  colnames(at) <- aa[,1]
  at <- at[order(rownames(at)), ,drop = FALSE]
  pdf(paste0(path, "summary/summary_heatmap.pdf"), width = 5, height = 6)
  p <- pheatmap(at, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 12)
  print(p)
  dev.off()
  
  # create summary boxplots
  pp <- melt(cbind(data.frame(sig_scores), type = type), id.var = "type")
  if(is.null(sig_meta)){
    labs <- as.character(unique(pp$variable))
    names(labs) <- labs
    pp$variable <- factor(as.character(pp$variable))
  } else {
    labs <- sig_meta
    pp$variable <- factor(pp$variable, names(sort(sig_meta)))
  }
  
  # plot dimensions
  if(n_sum_col == 2){
    wd <- 3.5
    ht <- 9
  } else {
    wd <- 5
    ht <- 6
  }
  
  pdf(paste0(path, "summary/summary_boxplot.pdf"), width = wd, height = ht)
  par(mar = c(2, 2, 2, 2))
  p <- ggplot(pp, aes(x=type, y=value, fill = type)) + 
    geom_boxplot(outlier.shape=NA, show.legend = FALSE) + 
    geom_jitter(color = "black", size = 0.2, alpha = 0.9, show.legend = FALSE) +
    facet_wrap(~ variable, scale="free_y", ncol=n_sum_col, labeller = as_labeller(labs), strip.position = "left") + scale_fill_manual(values = ann_col) + 
    my_theme() + 
    theme(strip.text.y = element_text(size = 12), text = element_text(size = 15), strip.background = element_blank(), strip.placement = "outside") +
    xlab("") + ylab("")
  print(p)
  dev.off()
}

#----------------------------------------------------------------------------------------
#' analyse_signature_by_trait
#' 
#' #' Generates scatter plots between a set of signature scores and scores for a single 
#' trait. Calcluates spearman p-values (excluding tied samples for either) trait
#' or signature, as well as Spearman's rho for the same samples. Corrects for multiple
#' testing using Holm's method
analyse_signature_by_trait <- function(trait, sig_scores, path, trait_id, alpha=0.05,
                                       sig_meta = NULL){
  plot_path <- paste0(path, "plots/")
  list_path <- paste0(path, "lists/")
  if(!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
  if(!dir.exists(list_path)) dir.create(list_path, recursive = TRUE)
  
  for(i in 1:ncol(sig_scores)){
    sig_id <- colnames(sig_scores)[i]
    message(sig_id)
    
    pdf(paste0(plot_path, sig_id, "_against_", trait_id, ".pdf"), width = 11, height = 9)
    p <- plot_scatter(trait = trait, y = sig_scores[,i], xlab = trait_id, ylab = sig_id)
    print(p)
    dev.off()
  }
  
  # Exclude tied samples
  trait_ties <- get_tie_ids(trait)
  sig_ties <- unlist(lapply(1:ncol(sig_scores), function(x) get_tie_ids(sig_scores[,x])))
  keep_samples <- !(1:length(trait) %in% union(sig_ties, trait_ties))
  p_vals <- sapply(1:ncol(sig_scores), function(x){
    cor.test(sig_scores[keep_samples,x], trait[keep_samples], method = "spearman")$p.value
  })
  rho_vals <- sapply(1:ncol(sig_scores), function(x){
    cor.test(sig_scores[keep_samples,x], trait[keep_samples], method = "spearman")$estimate
  })
  n_exc <- sum(!keep_samples)
  sig_tests <- data.frame(sig.id = colnames(sig_scores),
                          n.excluded.ties = n_exc,
                          n.total = length(trait),
                          rho = rho_vals,
                          p.val = p_vals, 
                          p.adj = p.adjust(p_vals, method = "holm"))
  sig_tests <- sig_tests[order(sig_tests$p.val), ]
  
  n_sig <- length(which(sig_tests$p.adj < 0.05))
  if(n_sig < 8){
    bs <- 8
  } else {
    bs <- 8
  }
  
  # plot summary barplot
  if(length(which(sig_tests$p.adj < alpha))){
    #pp <- subset(sig_tests, p.adj < 0.05)  
    pp <- sig_tests
    if(is.null(sig_meta)){
      pp$group <- factor(pp$sig.id, levels = rev(sort(pp$sig.id)))  
    } else {
      pp$group <- factor(sig_meta[as.character(pp$sig.id)], levels = rev(sort(sig_meta)))
    }
    
    # adjust plot size depending on number of rows
    pdf(paste0(path, trait_id, "_summary.pdf"), 1.5, 3, pointsize = 12)
    p <- ggplot(pp, aes(x = group, y = rho)) + geom_col(width=0.8) + 
      geom_text(aes(label = round(p.adj, 3), vjust = -0.2)) + 
      coord_flip() + my_theme2(base_size = bs) + ylim(-1,1) +
      ylab("Spearman's Rho")
    print(p)
    dev.off()
  }
  
  write.table(sig_tests, file = paste0(list_path, trait_id, "_correlation_test_p_values.csv"), sep = ",", col.names = NA)
}

#----------------------------------------------------------------------------------------
#' analyse_grouping_by_type
#' 
#' For each grouping, calculate cross-tab with types, and plot heatmap. Test for
#' significant deviation from expectation (under no assocation) between categories
#' using chi-squared test
analyse_grouping_by_type <- function(type, groupings, path, category_id, adjust_p_vals = TRUE, alpha = 0.05){
  plot_path <- paste0(path, "plots/")
  list_path <- paste0(path, "lists/")
  if(!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
  if(!dir.exists(list_path)) dir.create(list_path)
  
  for(i in 1:ncol(groupings)){
    grouping_id <- colnames(groupings)[i]
    message(grouping_id)
    
    cross_tab <- table(type, groupings[,i])
    
    pdf(paste0(plot_path, grouping_id, "_by_", category_id, ".pdf"))
    p <- pheatmap(apply(cross_tab, 2, function(x) x/sum(x)), cluster_rows = FALSE, cluster_cols = FALSE)
    print(p)
    dev.off()
  }
  
  p_vals <- apply(groupings, 2, function(x){
    # Use rule of thumb to check chi-squared test is appropriate
    stopifnot(validate_chisq(cross_tab)[[2]])
    cross_tab <- as.matrix(table(x, type))
    chisq.test(cross_tab)$p.value
  })
  
  grp_tests <- data.frame(p.val = p_vals)
  if(adjust_p_vals){
    grp_tests$p.adj = p.adjust(p_vals, method = "holm")
  } 
  grp_tests <- grp_tests[order(grp_tests$p.val), , drop = FALSE]
  
  write.table(grp_tests, file = paste0(list_path, category_id, "_chisq_test_p_values.csv"), sep = ",", col.names = NA)
  
  if(adjust_p_vals){
    signif_tests <- grp_tests[which(grp_tests$p.adj < alpha), ]  
  } else {
    signif_tests <- grp_tests[which(grp_tests$p.val < alpha), , drop=FALSE]
  }
  
  if(nrow(signif_tests) > 0){
    for(i in 1:nrow(signif_tests)){
      grp_id <- rownames(signif_tests)[i]
      grp_vals <- groupings[,match(grp_id, colnames(groupings))]
      cross_tab <- as.matrix(table(type, grp_vals))
      res <- chisq.posthoc.test(cross_tab, method = "bonferroni")
      write.table(res, file = paste0(list_path,  grp_id, "_by_", category_id, "_chisq_test_post_hoc.csv"),  sep = ",", col.names = NA)
    }  
  }
}

#----------------------------------------------------------------------------------------
#' analyse_survival
#' 
#' Wrapper survival analysis function. Generates Kaplan-Meier plot and runs associated log rank test for
#' every variable in var_ids, using quantized data for numeric variables. Runs CoxPH analyses for each
#' variable in var_ids with and without the covariates in cov_ids1 and cov_ids2. If there is more than 
#' one element in var_ids also runs combined Cox PH analyses using all var_ids elements at once, with 
#' and without covariates
#' Note: If cov_ids2 is non-empty, cov_ids1 must be as well
analyse_survival <- function(status_id, km_y_lab, survival_type_id, time_id, var_ids, f_surv_data, out_path, 
                             cov_ids1 = c(), cov_ids2 = c(), km_colors = NULL, km_x_lab = "Time (Months)", 
                             samples_id = "tnbc", fs = 3){
  stopifnot(is.data.frame(f_surv_data))
  stopifnot(length(cov_ids1) > 0 | length(cov_ids2) == 0)
  
  # create target paths
  km_path <- create_path(file.path(out_path, "kaplan_meier_plots/"))
  cox_path <- create_path(file.path(out_path, "coxph/"))
  if(length(cov_ids1) > 0) {
    cox_cov1_path <- create_path(file.path(out_path, "coxph_covariates/"))
  }
  if(length(cov_ids2) > 0){
    cox_cov2_path <- create_path(file.path(out_path, "coxph_extended_covariates/"))
  }
  
  # Kaplan-Meier analysis
  for(var_id in var_ids){
    quant_var_id <- c(var_id, paste0(var_id, ".quantile"))[is.numeric(f_surv_data[[var_id]]) + 1]
    srv <- drop_missing(f_surv_data[,c(time_id, status_id, quant_var_id)])
    km_fname <- paste0(km_path, var_id, "_", samples_id, "_", survival_type_id, "_km_survival.pdf")
    pdf(km_fname, width = 4, height = 4)
    p <- km_plot(status_id, km_y_lab, time_id, quant_var_id, srv, km_colors = km_colors, 
                 km_x_lab = km_x_lab, fs = fs)
    print(p)
    dev.off()
  }
  
  # separate CoxPH tests for var_ids without covariates
  for(var_id in var_ids){
    srv <- drop_missing(f_surv_data[, c(time_id, status_id, var_id)])
    cox_fname <- paste0(cox_path, var_id, "_", samples_id, "_", survival_type_id, "_coxph_res.csv")
    write.table(cox_regression(status_id, time_id, var_id, srv), cox_fname, sep = ",", col.names = NA)
    
    # separate CoxPH tests for var_ids with covariates
    if(length(cov_ids1) == 0) next
    srv <- drop_missing(f_surv_data[, c(time_id, status_id, var_id, cov_ids1)])
    cox_out <- cox_regression(status_id, time_id, var_id, srv, cov_ids = cov_ids1)
    cox_fname <- paste0(cox_cov1_path, var_id, "_", samples_id, "_", survival_type_id, "_cov_coxph_res.csv")
    write.table(cox_out, cox_fname, sep = ",", col.names = NA)
    
    # separate CoxPH tests for var_ids with extended covariates
    if(length(cov_ids2) == 0) next
    srv <- drop_missing(f_surv_data[, c(time_id, status_id, var_id, cov_ids2)])
    cox_out <- cox_regression(status_id, time_id, var_id, srv, cov_ids = cov_ids2)
    cox_fname <- paste0(cox_cov2_path, var_id, "_", samples_id, "_", survival_type_id, "_ext_cov_coxph_res.csv")
    write.table(cox_out, cox_fname, sep = ",", col.names = NA)
  }

  # do joint tests if there are multiple var_ids
  if(length(var_ids) > 1){
    # multivariate test without covariates
    srv <- drop_missing(f_surv_data[, c(time_id, status_id, var_ids)])
    cox_out <- cox_regression(status_id, time_id, var_ids, srv)
    cox_fname <- paste0(cox_path, "all_", samples_id, "_", survival_type_id, "_coxph_res.csv")
    write.table(cox_out, cox_fname, sep = ",", col.names = NA)
    
    # multivariate test with covariates
    if(length(cov_ids1) > 0) {
      srv <- drop_missing(f_surv_data[, c(time_id, status_id, var_ids, cov_ids1)])
      cox_out_cov <- cox_regression(status_id, time_id, var_ids, srv, cov_ids = cov_ids1)
      cox_fname_cov <- paste0(cox_cov1_path, "all_", samples_id, "_", survival_type_id, "_cov_coxph_res.csv")
      write.table(cox_out_cov, cox_fname_cov, sep = ",", col.names = NA)
    }
    
    # multivariate test with extended covariates
    if(length(cov_ids2) > 0) {
      srv <- drop_missing(f_surv_data[, c(time_id, status_id, var_ids, cov_ids2)])
      cox_out_cov <- cox_regression(status_id, time_id, var_ids, srv, cov_ids = cov_ids2)
      cox_fname_cov <- paste0(cox_cov2_path, "all_", samples_id, "_", survival_type_id, "_ext_cov_coxph_res.csv")
      write.table(cox_out_cov, cox_fname_cov, sep = ",", col.names = NA) 
    }
  }
}

#----------------------------------------------------------------------------------------
#' analyse_survival_multiple_testing
#' 
#' Wrapper survival analysis function. Generates Kaplan-Meier plot without log rank test for every 
#' variable in var_ids, using quantized data for numeric variables. Runs CoxPH analyses for each
#' variable in var_ids without. CoxPH p-values are corrected for multiple testing using Holm's method

#analyse_survival_multiple_testing("dss_status", sl["dss"], "dss", "dss_months", sids, 
#                                  surv_data, sp, cov_ids1 = cov_ids1, cov_ids2 = cov_ids2,
#                                  alpha = alpha)

analyse_survival_multiple_testing <- function(status_id, km_y_lab, survival_type_id, time_id, var_ids, f_surv_data,
                                              out_path, cov_ids1 = c(), cov_ids2 = c(), km_colors = NULL, 
                                              km_x_lab = "Time (Months)", samples_id = "tnbc", fs = 3, alpha = 0.05){
  stopifnot(is.data.frame(f_surv_data))
  stopifnot(length(cov_ids1) > 0 | length(cov_ids2) == 0)
  
  km_path <- create_path(file.path(out_path, "kaplan_meier_plots/"))
  cox_path <- create_path(file.path(out_path, "coxph/"))
  if(length(cov_ids1) > 0){
    cox_cov1_path <- create_path(file.path(out_path, "coxph_covariates/"))
  }
  if(length(cov_ids2) > 0){
    cox_cov2_path <- create_path(file.path(out_path, "coxph_extended_covariates/"))
  }
  
  # km plots
  for(var_id in var_ids){
    quant_var_id <- c(var_id, paste0(var_id, ".quantile"))[is.numeric(f_surv_data[[var_id]]) + 1]
    srv <- drop_missing(f_surv_data[,c(time_id, status_id, quant_var_id)])
    km_fname <- paste0(km_path, var_id, "_", samples_id, "_", survival_type_id, "_km_survival.pdf")
    pdf(km_fname, width = 4, height = 4)
    p <- km_plot(status_id, km_y_lab, time_id, quant_var_id, srv, km_colors = km_colors, 
                 km_x_lab = km_x_lab, fs = fs, p_val = FALSE, p_val_method = FALSE)
    print(p)
    dev.off()
  }
  
  # log-rank and Cox PH tests
  cox_list <- lapply(var_ids, function(x){
    srv <- drop_missing(f_surv_data[, c(time_id, status_id, x)])
    cox_regression(status_id, time_id, x, srv)
  })
  
  # output results
  cox_vals <- do.call(rbind, cox_list)
  cox_vals$p_adj <- p.adjust(cox_vals$`Pr(>|z|)`, method = "holm")
  p_res <- cbind(id = var_ids, cox_vals)
  p_res <- p_res[order(p_res$p_adj), ]
  
  res_fname <- paste0(cox_path, "all_", samples_id, "_", survival_type_id, "_coxph_res_holm.csv")
  write.table(p_res, res_fname, sep = ",", row.names = FALSE)
  
  sig_tests <- subset(p_res, p_adj < alpha)
  
  # test variables with significant associations with covariates included
  if(nrow(sig_tests) > 0 & length(cov_ids1) > 0){
    for(i in 1:nrow(sig_tests)){
      x <- sig_tests$id[i]
      srv <- drop_missing(f_surv_data[, c(time_id, status_id, x, cov_ids1)])
      cox_out <- cox_regression(status_id, time_id, x, srv, cov_ids = cov_ids1)
      cox_fname <- paste0(cox_cov1_path, x, "_", samples_id, "_", survival_type_id, "_cov_coxph_res.csv")
      write.table(cox_out, cox_fname, sep = ",", col.names = NA)
    }
  }
  
  # test variables with significant associations with extended covariates included
  if(nrow(sig_tests) > 0 & length(cov_ids2) > 0){
    for(i in 1:nrow(sig_tests)){
      x <- sig_tests$id[i]
      srv <- drop_missing(f_surv_data[, c(time_id, status_id, x, cov_ids2)])
      cox_out <- cox_regression(status_id, time_id, x, srv, cov_ids = cov_ids2)
      cox_fname <- paste0(cox_cov2_path, x, "_", samples_id, "_", survival_type_id, "_ext_cov_coxph_res.csv")
      write.table(cox_out, cox_fname, sep = ",", col.names = NA)
    }
  }
}

#----------------------------------------------------------------------------------------
#' analyse_primary_data_correlations
#' 
#' Primary function to look at correlations between covariates of interest in the data
analyse_primary_data_correlations <- function(exp_data, 
                                              exp_metadata,
                                              imm_scores,
                                              imm_metadata,
                                              imm_quantiles,
                                              patient_metadata,
                                              data_id,
                                              path_prefix = "analysis/human_primary_data/",
                                              alpha = 0.05,
                                              n_sum_col = 2,
                                              cov_ids1 = c(),
                                              cov_ids2 = c(),
                                              test_os = TRUE,
                                              test_dss = TRUE,
                                              test_pfs = FALSE){
  out_path <- paste0(path_prefix, data_id, "/")
  if(!dir.exists(out_path)) dir.create(out_path)
  
  # Analyse PRRX1 regulation
  prrx1_reg_plot_path <- paste0(out_path, "prrx1_regulation/plots/")
  prrx1_reg_list_path <- paste0(out_path, "prrx1_regulation/lists/")
  if(!dir.exists(prrx1_reg_plot_path)) dir.create(prrx1_reg_plot_path, recursive = TRUE)
  if(!dir.exists(prrx1_reg_list_path)) dir.create(prrx1_reg_list_path, recursive = TRUE)
  
  for(i in 1:length(prrx1_gsl)){
    message(names(prrx1_gsl)[i])
    assess_gsl_regulation(exp_data = exp_data, 
                          regulator = "PRRX1", 
                          gsl = prrx1_gsl[i], 
                          plot_path = prrx1_reg_plot_path, 
                          list_path = prrx1_reg_list_path)  
  }
  
  for(i in 1:length(prrx1_gml)){
    message(names(prrx1_gml)[i])
    assess_gml_regulation(exp_data = exp_data,
                          regulator = "PRRX1",
                          gml = prrx1_gml[i],
                          plot_path = prrx1_reg_plot_path,
                          list_path = prrx1_reg_list_path)
  }
  
  # Correlate prrx1 expression and targets with TNBC types
  prrx1_plot_path <- paste0(out_path, "prrx1_and_tnbc_type/plots/")
  prrx1_list_path <- paste0(out_path, "prrx1_and_tnbc_type/lists/")
  if(!dir.exists(prrx1_plot_path)) dir.create(prrx1_plot_path, recursive = TRUE)
  if(!dir.exists(prrx1_list_path)) dir.create(prrx1_list_path, recursive = TRUE)
  
  comp_ids <- c("prrx1_by_tnbc_type", 
                "prrx1_hs578_rna_targets_by_tnbc_type",
                "prrx1_mes_rna_targets_by_tnbc_type")
  y_lab_vals <- c("PRRX1 Expression",
                  "Hs578 RNA Targets",
                  "Mesenchymal RNA Targets")
  y_list <- list(as.numeric(exp_data["PRRX1",]),
                 exp_metadata$hs578_rna_targets,
                 exp_metadata$mes_rna_targets)
  
  for(i in 1:length(comp_ids)){
    message(comp_ids[i])
    
    pdf(paste0(prrx1_plot_path, comp_ids[i], ".pdf"), width = 11, height = 9)
    p <- plot_vioplot(type = exp_metadata$tnbc.type, y = y_list[[i]], ylab = y_lab_vals[i])
    print(p)
    dev.off()
    
    write_kruskal_test(type = exp_metadata$tnbc.type, y = y_list[[i]], path = prrx1_list_path, comp_id = comp_ids[i])
  }

  # Evaluate correlations between TNBC types Immune signatures
  analyse_signature_by_type(type = exp_metadata$tnbc.type, 
                            sig_scores = imm_scores, 
                            path = paste0(out_path, "immune_signatures_and_tnbc_types/"), 
                            category_id = "tnbc_type",
                            alpha = alpha,
                            sig_meta = imm_metadata,
                            n_sum_col = n_sum_col)
  
  # Evaluate correlations between PRRX1 expression and targets and Immune signatures
  analyse_signature_by_trait(trait = as.numeric(exp_data["PRRX1",]),
                             sig_scores = imm_scores,
                             path = paste0(out_path, "immune_signatures_and_prrx1/"),
                             trait_id = "prrx1_expression",
                             alpha = alpha,
                             sig_meta = imm_metadata)
  
  analyse_signature_by_trait(trait = exp_metadata$hs578_rna_targets,
                             sig_scores = imm_scores,
                             path = paste0(out_path, "immune_signatures_and_prrx1/"),
                             trait_id = "hs578_rna_targets",
                             alpha = alpha)
  
  analyse_signature_by_trait(trait = exp_metadata$mes_rna_targets,
                             sig_scores = imm_scores,
                             path = paste0(out_path, "immune_signatures_and_prrx1/"),
                             trait_id = "mes_rna_targets",
                             alpha = alpha)
  
  # Annotate survival data
  stopifnot(all(rownames(patient_metadata) == rownames(exp_metadata)))
  surv_data <- cbind(imm_scores, imm_quantiles, exp_metadata, patient_metadata)
  sl <- c(os = "Overall Survival Probability", 
          dss = "Disease-Specific Survival Probability",
          pfs = "Progression-Free Survival Probability")
  
  # Analyse survival associations with PRRX1 and targets
  surv_ids <- c("prrx1",
                "hs578_rna_targets",
                "mes_rna_targets")
  sp <- paste0(out_path, "prrx1_survival/")

  for(i in surv_ids){
    if(test_os)  analyse_survival("os_status",  sl["os"],  "os",  "os_months",  i, surv_data, sp,
                                  cov_ids1 = cov_ids1, cov_ids2 = cov_ids2)
    if(test_dss) analyse_survival("dss_status", sl["dss"], "dss", "dss_months", i, surv_data, sp,
                                  cov_ids1 = cov_ids1, cov_ids2 = cov_ids2)
    if(test_pfs) analyse_survival("pfs_status", sl["pfs"], "pfs", "pfs_months", i, surv_data, sp,
                                  cov_ids1 = cov_ids1, cov_ids2 = cov_ids2)
  }
  
  # Analyse survival associations with TNBC types
  sid <- "tnbc.type"
  sp <- paste0(out_path, "tnbc_type_survival/")
  kmc <- line_colours$type2
  
  if(test_os)  analyse_survival("os_status",  sl["os"],  "os",  "os_months",  sid,surv_data, sp, 
                                cov_ids1 = cov_ids1, cov_ids2 = cov_ids2, km_colors = kmc)
  if(test_dss) analyse_survival("dss_status", sl["dss"], "dss", "dss_months", sid, surv_data, sp, 
                                cov_ids1 = cov_ids1, cov_ids2 = cov_ids2, km_colors = kmc)
  if(test_pfs) analyse_survival("pfs_status", sl["pfs"], "pfs", "pfs_months", sid, surv_data, sp, 
                                cov_ids1 = cov_ids1, cov_ids2 = cov_ids2, km_colors = kmc)
  
  # Analyse survival assocations with immune signatures
  sids <- colnames(imm_scores)
  sp <- paste0(out_path, "immune_signatures_survival/")
  
  if(test_os)  analyse_survival_multiple_testing("os_status",  sl["os"],  "os",  "os_months",  sids, 
                                                 surv_data, sp, cov_ids1 = cov_ids1, cov_ids2 = cov_ids2,
                                                 alpha = alpha)
  if(test_dss) analyse_survival_multiple_testing("dss_status", sl["dss"], "dss", "dss_months", sids, 
                                                 surv_data, sp, cov_ids1 = cov_ids1, cov_ids2 = cov_ids2,
                                                 alpha = alpha)
  if(test_pfs) analyse_survival_multiple_testing("pfs_status", sl["pfs"], "pfs", "pfs_months", sids,
                                                 surv_data, sp, cov_ids1 = cov_ids1, cov_ids2 = cov_ids2,
                                                 alpha = alpha)
}
