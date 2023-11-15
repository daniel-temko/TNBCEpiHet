library(reshape2)
library(ComplexHeatmap)
library(VennDiagram)

top_pos_and_neg_genes <- function(weights_mat, factor_num, threshold = 0.2, n_genes = 200){
  pos_ids <- order(weights_mat[, factor_num], decreasing = TRUE)[1:n_genes]
  neg_ids <- order(weights_mat[, factor_num])[1:n_genes]
  return(list(pos = rownames(weights_mat)[pos_ids], 
              neg = rownames(weights_mat)[neg_ids]))
}

correlation_heatmap <- function(exp_mat, genes, fs = 15, col_func = NULL){
  if(is.null(col_func)){
    col_func <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "#EEEEEE", "red"))
  }
  p <- Heatmap(cor(t(exp_mat[genes, ])),
               col = col_func,
               cluster_rows = FALSE, cluster_columns = FALSE,
               show_row_names = FALSE,
               show_column_names = FALSE,
               heatmap_legend_param = list(title = "Pearson\nCorrelation",
                                           title_gp = gpar(fontsize = fs),
                                           labels_gp = gpar(fontsize = fs)))
  return(p)
}

pairwise_correlations <- function(exp_mat, genes){
  cor_mat <- cor(t(exp_mat[genes, , drop = FALSE]))
  return(cor_mat[upper.tri(cor_mat)])
}

cross_correlations <- function(exp_mat, genes1, genes2){
  stopifnot(length(genes1) == length(genes2))
  cor_mat <- cor(t(exp_mat[genes1, , drop = FALSE]), t(exp_mat[genes2, , drop = FALSE]))
  return(cor_mat[upper.tri(cor_mat)])
}

ecdf_vals <- function(vec, q_vals = seq(-1, 1, 0.01)){
  if(length(vec) == 0){
    return(NA)
  } else{
    return((ecdf(vec))(q_vals))  
  }
}

test_correlations_wrapper <- function(top_gene_list, exp_mat, n_rand_genes = 200,
                                      n_resample = 1000, n_threads = 20){
  require(foreach)
  require(doMC)
  registerDoMC(n_threads)
  stopifnot(nrow(exp_mat) > 2*n_rand_genes)
  
  cor_list <- lapply(top_gene_list, function(x){
    list(pos = pairwise_correlations(exp_mat, x$pos),
         neg = pairwise_correlations(exp_mat, x$neg),
         cross = cross_correlations(exp_mat, x$pos, x$neg))
  })
  
  mean_cor_vals <- sapply(cor_list, function(x){
    mean(c(x$pos, x$neg)) - mean(x$cross)
  })
  
  ecdf_list <- lapply(cor_list, function(x){
    list(pos = ecdf_vals(x$pos),
         neg = ecdf_vals(x$neg),
         all = ecdf_vals(c(x$pos, x$neg)),
         all_cross = ecdf_vals(x$cross))
  })
  
  rand_mean_cor_vals = unlist(foreach(i = 1:n_resample) %dopar% {
    set.seed(i)
    rand_sig <- sample(rownames(exp_mat), size = nrow(exp_mat), replace = FALSE)
    pos_genes <- rand_sig[1:n_rand_genes]
    neg_genes <- rev(rand_sig)[1:n_rand_genes]
    bg_cor_list <- list(pos = pairwise_correlations(exp_mat, pos_genes),
                        neg = pairwise_correlations(exp_mat, neg_genes),
                        cross = cross_correlations(exp_mat, pos_genes, neg_genes))
    mean(c(bg_cor_list$pos, bg_cor_list$neg)) -  mean(bg_cor_list$cross)
  })
  
  p_vals <- sapply(mean_cor_vals, function(x){
    length(which(rand_mean_cor_vals >= x)) / length(rand_mean_cor_vals)
  })
  
  set.seed(n_resample + 1)
  rand_sig <- sample(rownames(exp_mat), size = nrow(exp_mat), replace = FALSE)
  pos_genes <- rand_sig[1:n_rand_genes]
  neg_genes <- rev(rand_sig)[1:n_rand_genes]
  bg_cor <- list(pos = pairwise_correlations(exp_mat, pos_genes),
                 neg = pairwise_correlations(exp_mat, neg_genes),
                 cross = cross_correlations(exp_mat, pos_genes, neg_genes))
  bg_ecdf <- list(pos = ecdf_vals(bg_cor$pos),
                  neg = ecdf_vals(bg_cor$neg),
                  all = ecdf_vals(c(bg_cor$pos, bg_cor$neg)),
                  all_cross = ecdf_vals(bg_cor$cross))
  
  list(factors = list(mean_cor_vals = mean_cor_vals,
                      p_vals = p_vals,
                      ecdf_lists = ecdf_list),
       rand = list(mean_cor_vals = rand_mean_cor_vals,
                   ecdf = bg_ecdf))
}

plot_all_ecdfs <- function(f_test_res, pal){
  par(mar = c(6, 6, 6, 6))
  x_vals <- seq(-1, 1, 0.01)
  n_factors <- length(f_test_res$factors$ecdf_lists)
  plot(x_vals, f_test_res$rand$ecdf$all, col = "grey", type = "l", lwd = 2, lty = 2,
       xlab = "Pearson Correlation", ylab = expression(hat(F)(x)), cex.lab = 1.5,
       cex.axis = 1.5)
  for(i in 1:n_factors){
    par(new = TRUE)
    plot(x_vals, f_test_res$factors$ecdf_list[[i]]$all, col = pal[i], type = "l", lwd = 2, axes = FALSE,
         xlab = "", ylab = "")
  }
  legend("topleft", 
         legend = c("Random Genes", paste0("F", 1:n_factors, " Genes")), 
         col = c("grey", pal),
         lty = c(2, rep(1, n_factors)), lwd = 2,
         cex = 1.2)
  abline(v = 0, lty = 2)
}

plot_cross_ecdfs <- function(f_test_res, pal){
  par(mar = c(6, 6, 6, 6))
  x_vals <- seq(-1, 1, 0.01)
  n_factors <- length(f_test_res$factors$ecdf_lists)
  plot(x_vals, f_test_res$rand$ecdf$all_cross, col = "grey", type = "l", lwd = 2, lty = 2,
       xlab = "Pearson Correlation", ylab = expression(hat(F)(x)), cex.lab = 1.5,
       cex.axis = 1.5)
  for(i in 1:n_factors){
    par(new = TRUE)
    plot(x_vals, f_test_res$factors$ecdf_list[[i]]$all_cross, col = pal[i], type = "l", lwd = 2, axes = FALSE,
         xlab = "", ylab = "")
  }
  legend("bottomright", 
         legend = c("Random Genes", paste0("F", 1:n_factors, " Genes")), 
         col = c("grey", pal),
         lty = c(2, rep(1, n_factors)), lwd = 2,
         cex = 1.2)
  abline(v = 0, lty = 2)
}

#' Summarize results of permutation tests for all factors
plot_corr_significance <- function(f_test_res, pal){
  par(mar = c(6,6,6,6))
  p_val_df <- data.frame(mean_cor = f_test_res$factors$mean_cor_vals, p_val = f_test_res$factors$p_vals)
  p_val_df <- p_val_df[order(p_val_df$mean_cor),]
  x_vals <- c(test_res$rand$mean_cor_vals, test_res$factors$mean_cor_vals)
  x_min <- max(min(x_vals) - 0.05, -1)
  x_max <- min(max(x_vals) * 1.4, 1)
  n_factors <- length(f_test_res$factors$mean_cor_vals)
  n_resample <- length(f_test_res$rand$mean_cor_vals)
  x <- hist(f_test_res$rand$mean_cor_vals, col = "grey", border = "grey", 
            breaks = seq(x_min, x_max, 0.001), xlab = "Correlation Statistic",
            main = "", cex.lab = 1.5, cex.axis = 1.5)
  for(i in 1:n_factors){
    abline(v = f_test_res$factors$mean_cor_vals[i], lwd = 2, col = pal[i], lty = 2)
  }
  y_min <- min(x$counts)
  y_max <- max(x$counts)
  y_vals <- rev(seq(y_min, y_max, y_max/(n_factors)))[-(n_factors + 1)]
  # Add p-values
  for(i in 1:n_factors){
    x_val <- p_val_df$mean_cor[i]
    y_val <- y_vals[i]
    if(p_val_df$p_val[i] == 0){
      text(x_val, y_val, paste0("P<", 1/n_resample), cex = 1.2)
    } else {
      text(x_val, y_val, paste0("P=", p_val_df$p_val[i]), cex = 1.2)
    }
  }
  legend("topright", 
         legend = paste0("F", 1:n_factors), 
         col = pal, lwd = 2, cex = 1.2, lty = 2)
}

#' Create a 2-way Venn diagram for the overlap of the top_n features by absolute weight
#' in two MOFA models for a given view, for the facotrs specified by factor_new and factor_old,
#' in the new and old models. The weights for the old model are passed in weights_old 
#' and the weights for the new model are passed in weights_new
plot_overlap_venn <- function(weights_old, weights_new, factor_old, factor_new, 
                              view, fill_old, fill_new, top_n = 200, cex = 2){
  w_old <- weights_old[which((weights_old$view == view) & (weights_old$factor == factor_old)),]
  w_old <- w_old[order(abs(w_old$value), decreasing = TRUE), ]
  w_new <- weights_new[which((weights_new$view == view) & (weights_new$factor == factor_new)),]
  w_new <- w_new[order(abs(w_new$value), decreasing = TRUE), ]
  stopifnot(top_n < nrow(w_old))
  set_old <- w_old$feature[1:top_n]
  set_new <- w_new$feature[1:top_n]
  n_common <- length(intersect(set_old, set_new))
  draw.pairwise.venn(area1 = length(set_old), area2 = length(set_new), cross.area = n_common,
                     category = c(paste(factor_old, view), paste(factor_new, view)),
                     cat.pos = c(350, 10),
                     fill = c(fill_old, fill_new),
                     cat.cex = cex,
                     cex = cex)
}

#' Find the average ranking, by absolute weight of the common features that are in the top_n feautures
#' by absolute weight in two MOFA models, and return the results as a data.frame. Additionally, add
#' annotation on whether the feature was positively or negatively weighted in each of the models.
get_overlap_ranks <- function(weights_old, weights_new, factor_old, factor_new, view, top_n = 200){
  w_old <- weights_old[which((weights_old$view == view) & (weights_old$factor == factor_old)), ]
  w_old <- w_old[order(abs(w_old$value), decreasing = TRUE), ]
  w_new <- weights_new[which((weights_new$view == view) & (weights_new$factor == factor_new)), ]
  w_new <- w_new[order(abs(w_new$value), decreasing = TRUE), ]
  stopifnot(top_n < nrow(w_old))
  common_features <- intersect(w_old$feature[1:top_n], w_new$feature[1:top_n])
  rank_df <- data.frame(feature = common_features,
                        rank_old = match(common_features, w_old$feature),
                        rank_new = match(common_features, w_new$feature))
  rank_df$mean.rank <- (rank_df$rank_old + rank_df$rank_new) / 2
  rank_df$dir_old <- sapply(w_old$value[rank_df$rank_old], function(x){
    if(x > 0) {
      "pos"
    } else {
      "neg"
    }
  })
  rank_df$dir_new <- sapply(w_new$value[rank_df$rank_new], function(x){
    if(x > 0){
      "pos"
    } else {
      "neg"
    }
  })
  return(rank_df[order(rank_df$mean.rank, rank_df$feature), ])
}

#' Output a data.frame of p-values from hypergeometric tests, comparing the weights of each view
#' in each factor in an old MOFA model to the weights of every factor for the same view in a new 
#' MOFA model. The weights passed for each model are assumed to be filter to overlapping features
#' only in the calculation of the size of the feature universe for the hyper-geometric test. Each
#' test considers the overlap of the top_n weighted features by absolute weight for each model
compare_factor_weights_hypergeom <- function(weights_old, weights_new, top_n = 200){
  stopifnot(all(levels(weights_old$feature) == levels(weights_new$feature)))
  factors_old <- levels(weights_old$factor)
  factors_new <- levels(weights_new$factor)
  views <- levels(weights_old$view)
  comb_df <- expand.grid(view = views, new_factor = factors_new, old_factor = factors_old)
  comb_df <- comb_df[,c("old_factor", "new_factor", "view")]
  comb_df$p_val <- sapply(1:nrow(comb_df), function(x){
    old_ids <- which((weights_old$factor == comb_df$old_factor[x]) &
                       (weights_old$view == comb_df$view[x]))
    w_old <- droplevels(weights_old[old_ids, ])
    w_old <- w_old[order(abs(w_old$value), decreasing = TRUE), ]
    new_ids <- which((weights_new$factor == comb_df$new_factor[x]) & 
                       (weights_new$view == comb_df$view[x]))
    w_new <- droplevels(weights_new[new_ids, ])
    w_new <- w_new[order(abs(w_new$value), decreasing = TRUE), ]
    p_val <- hypergeom_p(w_old$feature, w_new$feature, top_n)
  })
  return(comb_df)
}

compare_factor_weights_pearson <- function(weights_old, weights_new){
  stopifnot(all(levels(weights_old$feature) == levels(weights_new$feature)))
  factors_old <- levels(weights_old$factor)
  factors_new <- levels(weights_new$factor)
  views <- levels(weights_old$view)
  comb_df <- expand.grid(view = views, new_factor = factors_new, old_factor = factors_old)
  comb_df <- comb_df[,c("old_factor", "new_factor", "view")]
  cor_list <- lapply(1:nrow(comb_df), function(x){
    old_ids <- which((weights_old$factor == comb_df$old_factor[x]) &
                       (weights_old$view == comb_df$view[x]))
    w_old <- weights_old[old_ids, ]
    new_ids <- which((weights_new$factor == comb_df$new_factor[x]) & 
                       (weights_new$view == comb_df$view[x]))
    w_new <- weights_new[new_ids, ]
    stopifnot(all(w_old$feature == w_new$feature))
    cor.test(w_old$value, w_new$value)
  })
  comb_df$cor <- sapply(cor_list, function(x) x$estimate)
  comb_df$p_val <- sapply(cor_list, function(x) x$p.value)
  return(comb_df)
}

#' Compute hyper-geometric test p-values for the overlap of the top_n
#' elements in vec1 and vec2
hypergeom_p <- function(vec1, vec2, top_n = 200){
  stopifnot(all(sort(vec1) == sort(vec2)))
  stopifnot(top_n < length(vec1))
  n_common <- length(intersect(vec1[1:top_n], vec2[1:top_n]))
  p_val <- phyper(n_common - 1, top_n, length(vec1) - top_n, top_n, lower.tail = FALSE)
  return(p_val)
}

#' Plot MOFA model weights, labeling the selected features given by label_features
plot_weights_3 <- function(mofa_object,
                           factor_num,
                           view,
                           label_features,
                           num_features = 10,
                           data_annotation = data_ann,
                           text_size = 20){
  #' Plot weights by rank and highlight features passed
  #' to label_features
  stopifnot(length(label_features) > 0)
  weights <- get_weights(mofa_object, 
                         factors = factor_num, 
                         views = view, 
                         as.data.frame = TRUE, 
                         scale = TRUE)
  feature_ids <- as.character(weights$feature)
  feature_names <- data_annotation[[view]][feature_ids, "mofa.name"]
  weights$feature_name <- feature_names
  weights$label <- ""
  lab_ids <- match(label_features, weights$feature)
  weights$label[lab_ids] <- weights$feature_name[lab_ids]
  weights <- weights[order(weights$value), ]
  weights$rank <- 1:nrow(weights)
  p <- ggplot(weights, aes(x = value, y = rank)) +
    geom_point() +
    geom_text_repel(aes(label = label),
                    max.overlaps = Inf,
                    force = 10,
                    direction = "both",
                    seed = 1,
                    ylim = c(nrow(weights) * 0.01, nrow(weights) * 0.99),
                    size = 6) +
    theme(text = element_text(size = text_size)) +
    xlab("Weight") +
    ylab("Rank") +
    xlim(c(-1, 1))
  return(p)
}