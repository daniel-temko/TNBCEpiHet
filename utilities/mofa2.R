source("../utilities/utils.R")
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(reshape2)
library(pheatmap)
library(DESeq2)
require(foreach)
library(doMC)
registerDoMC(22)

data_set_name_map <- c(methyl.se = "Methylation SE",
                       methyl.gb = "Methylation GB",
                       methyl.tss = "Methylation TSS",
                       chip.se = "ChIP Seq SE",
                       rna = "mRNA",
                       metab = "Metabolomics",
                       hms = "Histone MS")

#' Subset data to the top n_features features by variance, and sort in that order. If,
#' n_features is at least max_prop of the total features, then return the original matrix/
#' data.frame sorted by variance
SelectFeatures <- function(data, n_features, max_prop = 0.9){
  data_variances <- apply(data, 1, var)
  prop_features <- n_features / length(data_variances)
  if(prop_features >= max_prop){
    top_n <- length(data_variances)
  } else {
    top_n <- n_features
  }
  top_features <- names(sort(data_variances, decreasing = TRUE))[1:top_n]
  return(data[top_features, , drop = FALSE])
}

#' Reorders all columns of the matrices/data frames in data alphabetically and adds,
#' all NA columns for any sample absent in any entry. Converts elements to 
#' matrix
format_model_data <- function(data){
  data <- lapply(data, as.matrix)
  all_samples <- unique(unlist(lapply(data, colnames)))
  data <- lapply(data, function(x) {
    missing_cols <- all_samples[which(!all_samples %in% colnames(x))]
    if(length(missing_cols) > 0){
      na_mat <- matrix(NA, 
                       nrow = nrow(x), 
                       ncol = length(missing_cols), 
                       dimnames = list(NULL, missing_cols))
      cbind(x, na_mat)
    } else {
      x
    }
  })
  data <- lapply(data, function(x) x[,order(colnames(x))])
  return(data)
}

#' Load the mofa_objects given in model_metadata
load_mofa_objects <- function(model_metadata, data_path = "R_Data"){
  mofa_objects <- list()
  for(i in 1:length(model_metadata)){
    r_filename <- file.path(data_path, model_metadata[[i]]$r_file)
    mofa_objects[[names(model_metadata)[i]]] <- readRDS(r_filename)
  }
  return(mofa_objects)
}

#' Plot MOFA model residuals
plot_residuals <- function(mofa_object){
  data <- get_data(mofa_object)
  views <- names(data)
  factors <- get_factors(mofa_object)
  weights <- get_weights(mofa_object)
  
  data_exp <- lapply(views, function(x){
    sweep(weights[[x]] %*% t(factors$group1), 1, mofa_object@intercepts[[x]]$group1, "+")
  })
  
  plot_list <- lapply(1:length(data), function(x) {
    vals <- data[[x]]$group1
    keep <- apply(vals, 2, function(x) length(which(is.na(x))) == 0)
    vals <- vals[, keep]
    pred_vals <- data_exp[[x]][, keep]
    df <- data.frame(diff = as.numeric(vals - pred_vals))
    ggplot(df, aes(x = diff)) + 
      geom_histogram() +
      ggtitle(paste0(views[x])) +
      xlab("Residual")
  })
  plot_grid(plotlist = plot_list, ncol = 3)
}

#' Plot variance explained by each factor in each view
plot_variance_explained_raster <- function(mofa_object){
  df <- get_variance_explained(mofa_object, as.data.frame = TRUE)$r2_per_factor
  df$value <- df$value / 100
  n_factors <- mofa_object@dimensions$K
  factor_names <- paste0("Factor", 1:n_factors)
  plot_list <- lapply(factor_names, function(x){
    sub <- df[which(df$factor == x),]
    ggplot(sub, aes(x = view, y = value, fill = view)) +
      geom_bar(stat = "identity") +
      ylim(c(0, 0.5)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
  })
  plot_grid(plotlist = plot_list, ncol = 3)
}

# Plot the top weights for a given view and factor
plot_top_weights_2 <- function(object, 
                               view = 1, 
                               factor_num = 1, 
                               num_features = 10, 
                               abs = TRUE, 
                               scale = TRUE, 
                               sign = "all",
                               data_annotation = data_ann,
                               text_size = 20) 
{
  if (!is(object, "MOFA")) 
    stop("'object' has to be an instance of MOFA")
  if (num_features <= 0) 
    stop("'num_features' has to be greater than 0")
  if (sign == "all") {
    abs <- TRUE
  }
  if (is.numeric(view)) 
    view <- views_names(object)[view]
  stopifnot(view %in% views_names(object))
  view <- MOFA2:::.check_and_get_views(object, view)
  factors <- MOFA2:::.check_and_get_factors(object, factor_num)
  W <- get_weights(object, factors = factors, views = view, 
                   as.data.frame = TRUE)
  if (scale) 
    W$value <- W$value/max(abs(W$value))
  W <- W[W$value != 0, ]
  W$sign <- ifelse(W$value > 0, "+", "-")
  if (sign == "positive") {
    W <- W[W$value > 0, ]
  }
  else if (sign == "negative") {
    W <- W[W$value < 0, ]
  }
  if (abs) 
    W$value <- abs(W$value)
  W <- W[with(W, order(-abs(value))), ]
  W <- as.data.frame(top_n(group_by(W, factor), n = num_features, 
                           wt = value))
  feature_ids <- as.character(W$feature)
  feature_names <- data_annotation[[view]][feature_ids, "mofa.name"]
  W$feature_id <- feature_names
  if ((length(unique(W$view)) > 1) && (num_features > 0) && (any(duplicated(W[W$factor == 
                                                                           factors[1], ]$feature_id)))) {
    message("Duplicated feature names across views, we will add the view name as a prefix")
    W$feature_id <- paste(W$view, W$feature, sep = "_")
  }
  W$feature_id <- factor(W$feature_id, levels = rev(unique(W$feature_id)))
  p <- ggplot(W, aes_string(x = "feature_id", y = "value")) + 
    geom_point(size = 2) + geom_segment(aes_string(xend = "feature_id"), 
                                        size = 0.75, yend = 0) + scale_colour_gradient(low = "grey", 
                                                                                       high = "black") + coord_flip() + labs(y = "Weight") + 
    theme_bw() + theme(axis.title.x = element_text(color = "black"), 
                       axis.title.y = element_blank(), axis.text.y = element_text(size = rel(1.1), 
                                                                                  hjust = 1, color = "black"), axis.text.x = element_text(color = "black"), 
                       axis.ticks.y = element_blank(), axis.ticks.x = element_line(), 
                       legend.position = "top", legend.title = element_blank(), 
                       legend.text = element_text(color = "black"), legend.key = element_rect(fill = "transparent"), 
                       strip.text = element_text(size = rel(1.2)), panel.background = element_blank(), 
                       panel.spacing = unit(1, "lines"), panel.grid.major.y = element_blank(), 
    ) + facet_wrap(~factor, nrow = 1, scales = "free")
  if (sign == "negative") 
    p <- p + scale_x_discrete(position = "top")
  if (abs) {
    p <- p + ylim(0, max(W$value) + 0.1) + geom_text(label = W$sign, 
                                                     y = max(W$value) + 0.1, size = 10)
  }
  p <- p + theme(text = element_text(size = text_size))
  return(p)
}

#' Plot weights by rank and highlight num_features features with
#' top absolute weights
plot_weights_2 <- function(mofa_object,
                           factor_num,
                           view,
                           num_features = 10,
                           data_annotation = data_ann,
                           text_size = 20,
                           label = "top_absolute"){
  weights <- get_weights(mofa_object, 
                         factors = factor_num, 
                         views = view, 
                         as.data.frame = TRUE, 
                         scale = TRUE)
  feature_ids <- as.character(weights$feature)
  feature_names <- data_annotation[[view]][feature_ids, "mofa.name"]
  weights$feature_name <- feature_names
  weights$label <- ""
  if (label == "top_absolute"){
    lab_ids <- top_and_bottom_ids(abs(weights$value), 
                                  n_top = num_features,
                                  n_bottom = 0)$top_ids
  } else if (label == "top_pos_and_neg"){
    lab_ids <- unlist(top_and_bottom_ids(weights$value, 
                                         n_top = num_features, 
                                         n_bottom = num_features))
  }
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

#' Plot the samples ordinated by factors, and coloured according to type
plot_factor_raster <- function(mofa_object, 
                               factor_nums = 1:3, 
                               colour_by = "type", 
                               colour_map = line_colours){
  factors <- get_factors(mofa_object)$group1
  mofa_metadata <- samples_metadata(mofa_object)
  factor_names <- paste0("Factor", factor_nums)
  design <- expand.grid(factor_names, factor_names)
  types <- as.character(mofa_metadata[[colour_by]][match(rownames(factors), mofa_metadata$sample)])
  colours <- colour_map[[colour_by]][types]
  plot_list <- lapply(1:nrow(design), function(x){
    if(design$Var1[x] != design$Var2[x]){
      df <- data.frame(x = factors[, as.character(design$Var1[x])],
                       y = factors[, as.character(design$Var2[x])],
                       type = types)
      ggplot(df, aes(x = x, y = y, color = type)) + 
        geom_point(color = colours) +
        xlab(design$Var1[x]) +
        ylab(design$Var2[x])
    } else {
      NULL
    }
  })
  plot_grid(plotlist = plot_list, ncol = length(factor_nums))
}

#' Plot samples ordinated by factor scores for a single pair of factors, without 
#' labelling samples. Useful when there are many samples for readability
plot_unlabelled_factor_scatter <- function(mofa_object, factor_nums = c(2, 3), colour_by = "type", 
                                           colour_map = line_colours, text_size = 30){
  factors <- get_factors(mofa_object)$group1
  mofa_metadata <- samples_metadata(mofa_object)
  factor_names <- paste0("Factor", factor_nums)
  types <- as.character(mofa_metadata[[colour_by]][match(rownames(factors), mofa_metadata$sample)])
  colours <- colour_map[[colour_by]][types]
  df <- data.frame(x = factors[, as.character(factor_names[1])],
                   y = factors[, as.character(factor_names[2])],
                   type = types)
  ggplot(df, aes(x = x, y = y)) +
    geom_point(color = colours) +
    xlab(factor_names[1]) +
    ylab(factor_names[2]) +
    theme(text = element_text(size = text_size))
}

#' Plot samples ordinated by a single pair of factors, with sample labels
plot_labelled_factor_scatter <- function(mofa_object, factor_nums = c(2, 3), colour_by = "type", colour_map = line_colours){
  factors <- get_factors(mofa_object)$group1
  mofa_metadata <- samples_metadata(mofa_object)
  factor_names <- paste0("Factor", factor_nums)
  types <- as.character(mofa_metadata[[colour_by]][match(rownames(factors), mofa_metadata$sample)])
  colours <- colour_map[[colour_by]][types]
  df <- data.frame(sample = rownames(factors),
                   x = factors[, as.character(factor_names[1])],
                   y = factors[, as.character(factor_names[2])],
                   type = types)
  ggplot(df, aes(x = x, y = y)) +
    geom_point(color = colours) +
    xlab(factor_names[1]) +
    ylab(factor_names[2]) +
    geom_text_repel(label = df$sample,
                    max.overlaps = Inf)
}

#' Plot data heatmap showing the top feautres for a single factor and view.
#' Note: The sample_annotation argument is not used. Included for reverse compatability
plot_annotated_data_heatmap <- function(mofa_object, 
                                        factor_num, 
                                        view, 
                                        num_features = 50, 
                                        sample_annotation = line_annotation, 
                                        sample_colours = line_colours,
                                        data_annotation = data_ann,
                                        selection_method = "top_absolute"){
  mofa_metadata <- samples_metadata(mofa_object)
  factor_id <- paste0("Factor", factor_num)
  weights <- get_weights(mofa_object, factors = factor_id, views = view, as.data.frame = TRUE)
  factors <- get_factors(mofa_object, factors = factor_id, as.data.frame = TRUE)
  data <- get_data(mofa_object)[[view]]$group1
  weights <- weights[order(weights$value),]
  if (selection_method == "top_absolute"){
    select_ids <- top_and_bottom_ids(abs(weights$value), 
                                     n_top = num_features,
                                     n_bottom = 0)$top_ids
  } else if (selection_method == "top_pos_and_neg"){
    select_id_list <- top_and_bottom_ids(weights$value, 
                                         n_top = num_features, 
                                         n_bottom = num_features)
    select_ids <- c(select_id_list$top_ids, rev(select_id_list$bottom_ids))
  }
  weights <- weights[select_ids, ]
  
  # Remove NA columns
  keep <- apply(data, 2, function(x) length(which(is.na(x))) == 0)
  data <- data[, keep]
  
  # filter and rename features
  data <- data[as.character(weights$feature),]
  rownames(data) <- data_annotation[[view]][rownames(data), "mofa.name"]
  
  ann_col <- data.frame(row.names = colnames(data),
                        factor.value = factors$value[match(colnames(data), factors$sample)])
  sample_annotation <- mofa_metadata[rownames(ann_col), "type", drop = FALSE] 
  ann_col <- cbind(ann_col, sample_annotation)
  
  p <- pheatmap(data,
                annotation_col = ann_col,
                annotation_colors = sample_colours,
                cluster_rows = FALSE)
  return(p)
}

#' Plot raster of scatters between data totals and factor values
plot_sample_total_correlations <- function(mofa_object, 
                                           factor_num, 
                                           colour_by = "type", 
                                           colour_map = line_colours,
                                           name_map = data_set_name_map){
  factor_id <- paste0("Factor", factor_num)
  data <- get_data(mofa_object)
  mofa_metadata <- samples_metadata(mofa_object)
  factors <- get_factors(mofa_object, factors = factor_id, as.data.frame = TRUE)
  totals <- lapply(data, function(x){
    apply(x$group1, 2, sum)
  })
  plot_list <- lapply(1:length(totals), function(x){
    df <- data.frame(sample = factors$sample, factor.value = factors$value)
    df$total <- totals[[x]][df$sample]
    df <- na.omit(df)
    p <- ggplot(df, aes(x = factor.value, y = total)) + 
      ylab(paste0(name_map[names(totals)[x]], " Total")) +
      xlab(paste0("Factor ", factor_num))
    if(!is.na(colour_by)){
      df$type <- as.character(mofa_metadata[[colour_by]][match(df$sample, mofa_metadata$sample)])
      colours <- colour_map[[colour_by]][df$type]
      p <- p + geom_point(color = colours)
    } else {
      p <- p + geom_point()
    }
    p
  })
  plot_grid(plotlist = plot_list, ncol = 3)
}

#' Plot beeswarm for a single factor, with sample labeling and colouring
plot_labelled_factor_beeswarm <- function(mofa_object, 
                                          factor_num, 
                                          colour_by = "type", 
                                          colour_map = line_colours){
  factor_id <- paste0("Factor", factor_num)
  mofa_metadata <- samples_metadata(mofa_object)
  factors <- get_factors(mofa_object, factors = factor_id, as.data.frame = TRUE)
  df <- data.frame(sample = factors$sample, factor = factors$factor, factor.value = factors$value)
  df$type <- as.character(mofa_metadata[[colour_by]][match(df$sample, mofa_metadata$sample)])

  pos <- position_jitter(width = 0.3, seed = 20210105)
  p <- ggplot(df, aes(x = "factor", y = factor.value, color = type)) +
    geom_point(position = position_jitter(seed = 1)) +
    geom_text_repel(aes(label = sample), 
                    position = position_jitter(seed = 1),
                    max.overlaps = Inf,
                    size = 6) +
    scale_color_manual(values = colour_map[[colour_by]]) +
    theme(text = element_text(size = 30),
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ylab(paste0("Factor", factor_num)) +
    ggtitle(paste0("Factor", factor_num))
  return(p)
}

#' Output the GSEA results for a single factor at a given significance cut-off
#' in data.frame format
get_gsea_results <- function(gsea_object, factor_num, alpha = 0.05){
  factor_id <- paste0("Factor", factor_num)
  p_vals <- gsea_object$pval[, factor_id, drop = FALSE]
  p_adj <- gsea_object$pval.adj[, factor_id, drop = FALSE]
  stopifnot(all(rownames(p_vals) == rownames(p_adj)))
  out_df <- data.frame(set = rownames(p_vals), p.val = p_vals[,1], p.adj = p_adj[,1])
  out_df <- out_df[which(out_df$p.adj < alpha), ]
  out_df <- out_df[order(out_df$p.adj), ]
}

# Run GSEA on a mofa model using the MOFA GSEA function
run_gsea <- function(mofa_object, data_ann, factor_nums, view, gst_list, gsea_design){
  run_settings <- c("negative", "positive", "all")
  # Update mofa object names
  gene_col <- gsea_design[view, "gene.col.id"]
  feat_names <- as.character(features_names(mofa_object)[[view]])
  new_names <- as.character(data_ann[[view]][feat_names, gene_col])
  if(!all(feat_names == new_names)){
    features_names(mofa_object)[[view]] <- new_names
  }
  res <- list()
  for(i in 1:length(gst_list)){
    message(names(gst_list)[i])
    x <- list()
    for(j in 1:length(run_settings)){
      x[[j]] <- run_enrichment(mofa_object,
                     view = view,
                     factors = factor_nums,
                     feature.sets = gst_list[[i]],
                     sign = run_settings[j], 
                     statistical.test = "parametric")
    }
    names(x) <- c("neg", "pos", "all")
    res[[names(gst_list)[i]]] <- x
  }
  res
}

#' Wrapper function that analyses factors and weights in a MOFA model
#' and outputs summary plots to file
interrogate_mofa_model <- function(mofa_object, 
                                   data_ann, 
                                   model_id, 
                                   path_prefix = "analysis/mofa2/", 
                                   colour_list = NA, 
                                   colour_map=line_colours){
  out_path <- paste0(path_prefix, model_id, "/")
  dir.create(out_path)
  
  views <- names(mofa_object@data)
  num_factors <- mofa_object@dimensions$K
  
  # Plot model residuals as diagnostic
  pdf(paste0(out_path, "model_residuals.pdf"))
  p <- plot_residuals(mofa_object)
  print(p)
  dev.off()
  
  # Plot model overview
  pdf(paste0(out_path, "data_overview.pdf"))
  p <- plot_data_overview(mofa_object)
  print(p)
  dev.off()
  
  pdf(paste0(out_path, "factor_correlations.pdf"))
  p <- plot_factor_cor(mofa_object)
  print(p)
  dev.off()
  
  pdf(paste0(out_path, "variance_explained.pdf"))
  p <- plot_variance_explained(mofa_object, max_r2 = 15)
  print(p)
  dev.off()
  
  pdf(paste0(out_path, "total_variance_explained.pdf"))
  p <- plot_variance_explained(mofa_object, plot_total = TRUE, max_r2 = 15)[[2]]
  print(p)
  dev.off()
  
  pdf(paste0(out_path, "variance_explained_raster.pdf"))
  p <- plot_variance_explained_raster(mofa_object)
  print(p)
  dev.off()
  
  pdf(paste0(out_path, "factors1to3_raster.pdf"))
  p <- plot_factor_raster(mofa_object, colour_map = colour_map)
  print(p)
  dev.off()
  
  mofa_metadata <- samples_metadata(mofa_object)
  phenotype_covariates <- colnames(mofa_metadata)[grep("type|total", colnames(mofa_metadata))]
  tech_covariates <- colnames(mofa_metadata)[grep("batch", colnames(mofa_metadata))]
  
  pdf(paste0(out_path, "factor_phenotype_correlations.pdf"))
  p <- try(correlate_factors_with_covariates(mofa_object, covariates = phenotype_covariates), silent = TRUE)
  print(p)
  dev.off()
  
  pdf(paste0(out_path, "factor_technical_correlations.pdf"))
  p <- try(correlate_factors_with_covariates(mofa_object, covariates = tech_covariates), silent = TRUE)
  print(p)
  dev.off()
  
  # Interrogate factors
  for(i in 1:num_factors){
    factor_name <- paste0("factor", i)
    feature_plot_path <- paste0(out_path, factor_name, "/feature_plots/")
    sample_plot_path <- paste0(out_path, factor_name, "/sample_plots/")
    feature_list_path <- paste0(out_path, factor_name, "/feature_lists/")
    
    dir.create(feature_plot_path, recursive = TRUE)
    dir.create(feature_list_path, recursive = TRUE)
    dir.create(sample_plot_path, recursive = TRUE)
    
    # Sample plots
    pdf(paste0(sample_plot_path, factor_name, "_sample_total_correlations.pdf"))
    p <- plot_sample_total_correlations(mofa_object, factor_num = i, colour_by = NA, colour_map = colour_map)
    print(p)
    dev.off()
    
    pdf(paste0(sample_plot_path, factor_name, "_sample_total_correlations_by_type.pdf"))
    p <- plot_sample_total_correlations(mofa_object, factor_num = i, colour_map = colour_map)
    print(p)
    dev.off()
    
    pdf(paste0(sample_plot_path, factor_name, "_beeswarm.pdf"))
    p <- plot_labelled_factor_beeswarm(mofa_object, factor = i, colour_map = colour_map)
    print(p)
    dev.off()
    
    # Batch-based plots
    if("rna.batch" %in% colnames(mofa_metadata)){
      pdf(paste0(sample_plot_path, factor_name, "_beeswarm_by_rna_batch.pdf"))
      p <- plot_labelled_factor_beeswarm(mofa_object, factor = i, colour_by = "rna.batch", colour_map = colour_list)
      print(p)
      dev.off()  
    }
    
    if("methyl.batch" %in% colnames(mofa_metadata)){
      pdf(paste0(sample_plot_path, factor_name, "_beeswarm_by_methyl_batch.pdf"))
      p <- plot_labelled_factor_beeswarm(mofa_object, factor = i, colour_by = "methyl.batch", colour_map = colour_list)
      print(p)
      dev.off()  
    }
    
    if("metab.batch" %in% colnames(mofa_metadata)){
      pdf(paste0(sample_plot_path, factor_name, "_beeswarm_by_metab_batch.pdf"))
      p <- plot_labelled_factor_beeswarm(mofa_object, factor = i, colour_by = "metab.batch", colour_map = colour_list)
      print(p)
      dev.off()  
    }
    
    for(j in 1:length(views)){
      view <- views[j]
      view_name <- gsub("\\.", "_", view)
      message("Analysing ", view_name, " ", factor_name)
      
      # Feature plots
      pdf(paste0(feature_plot_path, "weights_top_abs_weights_highlighted_", factor_name, "_", view_name, ".pdf"))
      p <- plot_weights_2(mofa_object, factor_num = i, view = views[j], data_annotation = data_ann)
      print(p)
      dev.off()
      
      pdf(paste0(feature_plot_path, "weights_top_pos_and_neg_weights_highlighted_", factor_name, "_", view_name, ".pdf"))
      p <- plot_weights_2(mofa_object, factor_num = i, view = views[j], num_features = 5, label = "top_pos_and_neg", data_annotation = data_ann)
      print(p)
      dev.off()
      
      pdf(paste0(feature_plot_path, "top_abs_weights_", factor_name, "_", view_name, ".pdf"))
      p <- plot_top_weights_2(mofa_object, factor_num = i, view = views[j], data_annotation = data_ann)
      print(p)
      dev.off()
      
      pdf(paste0(feature_plot_path, "top_pos_weights_", factor_name, "_", view_name, ".pdf"))
      p <- plot_top_weights_2(mofa_object, factor_num = i, view = views[j], sign = "positive", data_annotation = data_ann)
      print(p)
      dev.off()
      
      pdf(paste0(feature_plot_path, "top_neg_weights_", factor_name, "_", view_name, ".pdf"))
      p <- plot_top_weights_2(mofa_object, factor_num = i, view = views[j], sign = "negative", data_annotation = data_ann)
      print(p)
      dev.off()
      
      pdf(file = paste0(feature_plot_path, "heatmap_top_abs_weights_", factor_name, "_", view_name, ".pdf"), width = 11, height = 9)
      p <- plot_annotated_data_heatmap(mofa_object, factor_num = i, view = views[j], sample_colours = colour_map, data_annotation = data_ann)
      print(p)
      dev.off()
      
      pdf(file = paste0(feature_plot_path, "heatmap_top_pos_and_neg_weights_", factor_name, "_", view_name, ".pdf"), width = 11, height = 9)
      p <- plot_annotated_data_heatmap(mofa_object, 
                                       factor_num = i, 
                                       view = views[j], 
                                       num_features = 25, 
                                       selection_method = "top_pos_and_neg",
                                       sample_colours = colour_map,
                                       data_annotation = data_ann)
      dev.off()
      
      # Write weights to file
      weights <- get_weights(mofa_object, factors = i, views = views[j], scale = TRUE, as.data.frame = TRUE)
      weights <- weights[order(abs(weights$value), decreasing = TRUE), ]
      feat_ann <- data_ann[[view]][as.character(weights$feature), ]
      weights <- cbind(weights[,c("feature", "value")], feat_ann)
      write.table(weights, paste0(feature_list_path, factor_name, "_", view_name, "_scaled_weights.csv"), sep = ",", row.names = FALSE)
    }
  }
}

#' Wrapper function to run GSEA analysis on a trained MOFA model
run_gsea_analysis <- function(mofa_object, 
                              data_ann, 
                              model_id, 
                              gst_list, 
                              gsea_design,
                              path_prefix = "analysis/mofa2/"){
  out_paths <- lapply(1:length(gst_list), function(x){
    paste0(path_prefix, model_id, "/", "gsea/", names(gst_list)[x], "/")  
  })
  for(i in 1:length(gst_list)){
    dir.create(out_paths[[i]], recursive = TRUE)
  }

  views <- names(mofa_object@data)
  num_factors <- mofa_object@dimensions$K
  
  # Run GSEA
  res <- foreach (i = 1:nrow(gsea_design)) %dopar% {
    run_gsea(mofa_object = mofa_object,
             data_ann = data_ann,
             factor_nums = 1:num_factors,
             view = rownames(gsea_design)[i],
             gst_list = gst_list,
             gsea_design = gsea_design)
  }
  names(res) <- rownames(gsea_design)
  
  for(i in 1:num_factors){
    factor_name <- paste0("factor", i)
    gene_set_list_paths <- lapply(1:length(gst_list), function(x) {
      paste0(out_paths[[x]], factor_name, "/")
    })
    for(j in 1:length(gst_list)){
      dir.create(gene_set_list_paths[[j]], recursive = TRUE)  
    }
    for(j in 1:nrow(gsea_design)){
      view <- rownames(gsea_design)[j]
      view_name <- gsub("\\.", "_", view)
      message(view_name, " ", factor_name)
      
      gsea_res <- res[[view]]
      
      # Write results
      for(k in 1:length(gst_list)){
        gst_id <- names(gsea_res)[k]
        out_df <- get_gsea_results(gsea_res[[gst_id]]$neg, factor_num = i)
        fn <- paste0(gene_set_list_paths[[k]], gst_id, "_neg_", factor_name, "_", view_name, ".csv")
        write.table(out_df, fn, sep = ",", row.names = FALSE)
        
        out_df <- get_gsea_results(gsea_res[[gst_id]]$pos, factor_num = i)
        fn <- paste0(gene_set_list_paths[[k]], gst_id, "_pos_", factor_name, "_", view_name, ".csv")
        write.table(out_df, fn, sep = ",", row.names = FALSE)
        
        out_df <- get_gsea_results(gsea_res[[gst_id]]$all, factor_num = i)
        fn <- paste0(gene_set_list_paths[[k]], gst_id, "_all_", factor_name, "_", view_name, ".csv")
        write.table(out_df, fn, sep = ",", row.names = FALSE)
      }
    }
  }
}
