library(sctransform)

## Adjusted version of Seurat SCTransform wrapper for imbalanced batches
## Behaviour should match default SCTransform in other respects
sct_adj_rescale <- function(sc, cor_df, restrict_cells, seedval = 1448145, 
                            variable.features.n=3000, verbose = TRUE){
  set.seed(seedval)
  
  sc.mat <- GetAssayData(object = sc, assay = "RNA", slot = "counts")
  
  clip.range <- c(-sqrt(ncol(sc.mat) / 30), sqrt(ncol(sc.mat) / 30))
  
  cell_attr <- sc@meta.data[,c("batch"), drop=FALSE]
  cell_attr$cell.line <- sc@meta.data$orig.ident
  cell_attr$umi <- sc@meta.data$nCount_RNA
  cell_attr$log_umi <- log10(cell_attr$umi)
  cell_attr <- cbind(cell_attr, cor_df)
  
  ## Fit model based on reference cell lines only
  vst_out <- vst_rescale(sc.mat, 
                         latent_var = c("log_umi"),
                         latent_var_nonreg = colnames(cor_df),
                         cell_attr = cell_attr,
                         return_gene_attr = TRUE, 
                         return_cell_attr = FALSE,
                         return_corrected_umi = FALSE,
                         verbosity = 2,
                         restrict_cells = restrict_cells)
  
  ## Correct counts for all cell lines
  sct_cnts <- correct_counts(x = vst_out, umi = sc.mat, cell_attr = cell_attr)
  sct_res <- get_residuals(vst_out = vst_out, umi = sc.mat, cell_attr = cell_attr)
  
  feature.variance <- get_residual_var(vst_out = vst_out, umi = sc.mat, cell_attr = cell_attr)
  feature.variance <- sort(x = feature.variance, decreasing = TRUE)
  top.features <- names(x = feature.variance)[1:min(variable.features.n, length(x = feature.variance))]
  
  sct_scale <- sct_res[top.features, ]
  sct_scale[sct_scale < clip.range[1]] <- clip.range[1]
  sct_scale[sct_scale > clip.range[2]] <- clip.range[2]
  sct_scale <- ScaleData(sct_scale, do.scale = FALSE, scale.max = Inf)
  
  sct_assay <- CreateAssayObject(counts = sct_cnts)
  VariableFeatures(object = sct_assay) <- top.features
  sct_assay <- SetAssayData(object = sct_assay, slot = "data", 
                            new.data = log1p(x = GetAssayData(object = sct_assay, 
                                                              slot = "counts")))
  sct_assay <- SetAssayData(sct_assay, slot = "scale.data", new.data = sct_scale)
  vst_out$y <- NULL
  vst_out$arguments$sct.clip.range <- clip.range
  Misc(object = sct_assay, slot = "vst.out") <- vst_out
  sct_assay[[paste0("asct.", names(x = vst_out$gene_attr))]] <- vst_out$gene_attr
  sct_assay[["asct.variable"]] <- rownames(x = sct_assay[[]]) %in% top.features
  sc[["asct"]] <- sct_assay
  if(verbose){
    message("Set default assay to asct")
  }
  DefaultAssay(sc) <- "asct"
  return(sc)
}

## Adjusted version of Seurat SCTransform wrapper for imbalanced batches
## Behaviour should match default SCTransform in other respects
sct_adj_rescale_pos <- function(sc, cor_df, restrict_cells, seedval = 1448145, 
                            variable.features.n=3000, verbose = TRUE){
  set.seed(seedval)
  
  sc.mat <- GetAssayData(object = sc, assay = "RNA", slot = "counts")
  
  clip.range <- c(-sqrt(ncol(sc.mat) / 30), sqrt(ncol(sc.mat) / 30))
  
  cell_attr <- sc@meta.data[,c("batch"), drop=FALSE]
  cell_attr$cell.line <- sc@meta.data$orig.ident
  cell_attr$umi <- sc@meta.data$nCount_RNA
  cell_attr$log_umi <- log10(cell_attr$umi)
  cell_attr <- cbind(cell_attr, cor_df)
  
  ## Fit model based on reference cell lines only
  vst_out <- vst_rescale(sc.mat, 
                         latent_var = c("log_umi"),
                         latent_var_nonreg = colnames(cor_df),
                         cell_attr = cell_attr,
                         return_gene_attr = TRUE, 
                         return_cell_attr = FALSE,
                         return_corrected_umi = FALSE,
                         verbosity = 2,
                         restrict_cells = restrict_cells,
                         pos_coefs_only = TRUE)
  
  ## Correct counts for all cell lines
  sct_cnts <- correct_counts(x = vst_out, umi = sc.mat, cell_attr = cell_attr)
  sct_res <- get_residuals(vst_out = vst_out, umi = sc.mat, cell_attr = cell_attr)
  
  feature.variance <- get_residual_var(vst_out = vst_out, umi = sc.mat, cell_attr = cell_attr)
  feature.variance <- sort(x = feature.variance, decreasing = TRUE)
  top.features <- names(x = feature.variance)[1:min(variable.features.n, length(x = feature.variance))]
  
  sct_scale <- sct_res[top.features, ]
  sct_scale[sct_scale < clip.range[1]] <- clip.range[1]
  sct_scale[sct_scale > clip.range[2]] <- clip.range[2]
  sct_scale <- ScaleData(sct_scale, do.scale = FALSE, scale.max = Inf)
  
  sct_assay <- CreateAssayObject(counts = sct_cnts)
  VariableFeatures(object = sct_assay) <- top.features
  sct_assay <- SetAssayData(object = sct_assay, slot = "data", 
                            new.data = log1p(x = GetAssayData(object = sct_assay, 
                                                              slot = "counts")))
  sct_assay <- SetAssayData(sct_assay, slot = "scale.data", new.data = sct_scale)
  vst_out$y <- NULL
  vst_out$arguments$sct.clip.range <- clip.range
  Misc(object = sct_assay, slot = "vst.out") <- vst_out
  sct_assay[[paste0("asct.", names(x = vst_out$gene_attr))]] <- vst_out$gene_attr
  sct_assay[["asct.variable"]] <- rownames(x = sct_assay[[]]) %in% top.features
  sc[["asct"]] <- sct_assay
  if(verbose){
    message("Set default assay to asct")
  }
  DefaultAssay(sc) <- "asct"
  return(sc)
}

vst_rescale <- function (umi, cell_attr = NULL, latent_var = c("log_umi"), batch_var = NULL, 
          latent_var_nonreg = NULL, n_genes = 2000, n_cells = NULL, 
          method = "poisson", do_regularize = TRUE, theta_regularization = "od_factor", 
          res_clip_range = c(-sqrt(ncol(umi)), sqrt(ncol(umi))), bin_size = 500, 
          min_cells = 5, residual_type = "pearson", return_cell_attr = FALSE, 
          return_gene_attr = TRUE, return_corrected_umi = FALSE, min_variance = -Inf, 
          bw_adjust = 3, gmean_eps = 1, theta_estimation_fun = "theta.ml", 
          theta_given = NULL, verbosity = 2, verbose = NULL, show_progress = NULL,
          restrict_cells = NULL, pos_coefs_only = FALSE) 
{
  arguments <- as.list(environment())
  arguments <- arguments[!names(arguments) %in% c("umi", "cell_attr")]
  if (!is.null(verbose)) {
    warning("The 'verbose' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)", 
            immediate. = TRUE, call. = FALSE)
    verbosity <- as.numeric(verbose)
  }
  if (!is.null(show_progress)) {
    warning("The 'show_progress' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)", 
            immediate. = TRUE, call. = FALSE)
    if (show_progress) {
      verbosity <- 2
    }
    else {
      verbosity <- min(verbosity, 1)
    }
  }
  if (method == "glmGamPoi") {
    glmGamPoi_check <- requireNamespace("glmGamPoi", quietly = TRUE)
    if (!glmGamPoi_check) {
      stop("Please install the glmGamPoi package. See https://github.com/const-ae/glmGamPoi for details.")
    }
  }
  if (startsWith(x = method, prefix = "offset")) {
    cell_attr <- NULL
    latent_var <- c("log_umi")
    batch_var <- NULL
    latent_var_nonreg <- NULL
    n_genes <- NULL
    n_cells <- NULL
    do_regularize <- FALSE
    if (is.null(theta_given)) {
      theta_given <- 100
    }
    else {
      theta_given <- theta_given[1]
    }
  }
  times <- list(start_time = Sys.time())
  cell_attr <- sctransform:::make_cell_attr(umi, cell_attr, latent_var, batch_var, 
                              latent_var_nonreg, verbosity)
  if (!is.null(batch_var)) {
    cell_attr[, batch_var] <- as.factor(cell_attr[, batch_var])
    batch_levels <- levels(cell_attr[, batch_var])
  }
  genes_cell_count <- rowSums(umi >= 0.01)
  genes <- rownames(umi)[genes_cell_count >= min_cells]
  umi <- umi[genes, ]
  genes_log_gmean <- log10(sctransform:::row_gmean(umi, eps = gmean_eps))
  if (!do_regularize && !is.null(n_genes)) {
    if (verbosity > 0) {
      message("do_regularize is set to FALSE, will use all genes")
    }
    n_genes <- NULL
  }
  if (!is.null(n_cells) && n_cells < ncol(umi)) {
    cells_step1 <- sample(x = colnames(umi), size = n_cells)
    if (!is.null(batch_var)) {
      dropped_batch_levels <- setdiff(batch_levels, levels(droplevels(cell_attr[cells_step1, 
                                                                                batch_var])))
      if (length(dropped_batch_levels) > 0) {
        stop("Dropped batch levels ", dropped_batch_levels, 
             ", set n_cells higher")
      }
    }
    genes_cell_count_step1 <- rowSums(umi[, cells_step1] > 
                                        0)
    genes_step1 <- rownames(umi)[genes_cell_count_step1 >= 
                                   min_cells]
    genes_log_gmean_step1 <- log10(row_gmean(umi[genes_step1, 
                                                 cells_step1], eps = gmean_eps))
  }
  else {
    cells_step1 <- colnames(umi)
    genes_step1 <- genes
    genes_log_gmean_step1 <- genes_log_gmean
  }
  data_step1 <- cell_attr[cells_step1, , drop = FALSE]
  if (!is.null(n_genes) && n_genes < length(genes_step1)) {
    log_gmean_dens <- density(x = genes_log_gmean_step1, 
                              bw = "nrd", adjust = 1)
    sampling_prob <- 1/(approx(x = log_gmean_dens$x, y = log_gmean_dens$y, 
                               xout = genes_log_gmean_step1)$y + .Machine$double.eps)
    genes_step1 <- sample(x = genes_step1, size = n_genes, 
                          prob = sampling_prob)
    genes_log_gmean_step1 <- log10(sctransform:::row_gmean(umi[genes_step1, 
                                                 cells_step1], eps = gmean_eps))
  }
  if (!is.null(batch_var)) {
    model_str <- paste0("y ~ (", paste(latent_var, collapse = " + "), 
                        ") : ", batch_var, " + ", batch_var, " + 0")
  }
  else {
    model_str <- paste0("y ~ ", paste(latent_var, collapse = " + "))
  }
  bin_ind <- ceiling(x = 1:length(x = genes_step1)/bin_size)
  max_bin <- max(bin_ind)
  if (verbosity > 0) {
    message("Variance stabilizing transformation of count matrix of size ", 
            nrow(umi), " by ", ncol(umi))
    message("Model formula is ", model_str)
  }
  times$get_model_pars = Sys.time()
  model_pars <- sctransform:::get_model_pars(genes_step1, bin_size, umi, 
                               model_str, cells_step1, method, data_step1, theta_given, 
                               theta_estimation_fun, verbosity)
  min_theta <- 1e-07
  if (any(model_pars[, "theta"] < min_theta)) {
    if (verbosity > 0) {
      msg <- sprintf("There are %d estimated thetas smaller than %g - will be set to %g", 
                     sum(model_pars[, "theta"] < min_theta), min_theta, 
                     min_theta)
      message(msg)
    }
    model_pars[, "theta"] <- pmax(model_pars[, "theta"], 
                                  min_theta)
  }
  times$reg_model_pars = Sys.time()
  if (do_regularize) {
    model_pars_fit <- sctransform:::reg_model_pars(model_pars, genes_log_gmean_step1, 
                                     genes_log_gmean, cell_attr, batch_var, cells_step1, 
                                     genes_step1, umi, bw_adjust, gmean_eps, theta_regularization, 
                                     verbosity)
    model_pars_outliers <- attr(model_pars_fit, "outliers")
  }
  else {
    model_pars_fit <- model_pars
    model_pars_outliers <- rep(FALSE, nrow(model_pars))
  }
  regressor_data <- model.matrix(as.formula(gsub("^y", "", 
                                                 model_str)), cell_attr)
  if (!is.null(latent_var_nonreg)) {
    if (verbosity > 0) {
      message("Estimating parameters for following non-regularized variables: ", 
              latent_var_nonreg)
    }
    if (!is.null(batch_var)) {
      model_str_nonreg <- paste0("y ~ (", paste(latent_var_nonreg, 
                                                collapse = " + "), ") : ", batch_var, " + ", 
                                 batch_var, " + 0")
    }
    else {
      model_str_nonreg <- paste0("y ~ ", paste(latent_var_nonreg, 
                                               collapse = " + "), " + 0")
    }
    if(is.null(restrict_cells)){
      cells_non_reg <- colnames(umi)
    } else{
      stopifnot(is.numeric(restrict_cells))
      cells_non_reg <- colnames(umi)[restrict_cells]
    }
    times$get_model_pars_nonreg = Sys.time()
    model_pars_nonreg <- get_model_pars_nonreg_rescale(genes, bin_size, 
                                      model_pars_fit, regressor_data[cells_non_reg, ], umi[, cells_non_reg], 
                                      model_str_nonreg, cell_attr[cells_non_reg, ], verbosity,
                                      latent_var_nonreg)
    ## adjust coeffs
    if(pos_coefs_only){
      tot_effect <- model_pars_nonreg[,2] - model_pars_nonreg[,1]
      model_pars_nonreg[,1] <- abs(tot_effect) * (sign(tot_effect) == -1)
      model_pars_nonreg[,2] <- abs(tot_effect) * (sign(tot_effect) ==  1)
    } else {
      model_pars_nonreg <- sweep(model_pars_nonreg, 1, rowMeans(model_pars_nonreg))  
    }
    
    regressor_data_nonreg <- model.matrix(as.formula(gsub("^y", 
                                                          "", model_str_nonreg)), cell_attr)
    model_pars_final <- cbind(model_pars_fit, model_pars_nonreg)
    regressor_data_final <- cbind(regressor_data, regressor_data_nonreg)
  }
  else {
    model_str_nonreg <- ""
    model_pars_nonreg <- c()
    model_pars_final <- model_pars_fit
    regressor_data_final <- regressor_data
  }
  times$get_residuals = Sys.time()
  if (!residual_type == "none") {
    if (verbosity > 0) {
      message("Second step: Get residuals using fitted parameters for ", 
              length(x = genes), " genes")
    }
    bin_ind <- ceiling(x = 1:length(x = genes)/bin_size)
    max_bin <- max(bin_ind)
    if (verbosity > 1) {
      pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
    }
    res <- matrix(NA_real_, length(genes), nrow(regressor_data_final), 
                  dimnames = list(genes, rownames(regressor_data_final)))
    for (i in 1:max_bin) {
      genes_bin <- genes[bin_ind == i]
      mu <- exp(tcrossprod(model_pars_final[genes_bin, 
                                            -1, drop = FALSE], regressor_data_final))
      y <- as.matrix(umi[genes_bin, , drop = FALSE])
      res[genes_bin, ] <- switch(residual_type, pearson = sctransform:::pearson_residual(y, 
                                                                           mu, model_pars_final[genes_bin, "theta"], min_var = min_variance), 
                                 deviance = sctransform:::deviance_residual(y, mu, model_pars_final[genes_bin, 
                                                                                      "theta"]), stop("residual_type ", residual_type, 
                                                                                                      " unknown - only pearson and deviance supported at the moment"))
      if (verbosity > 1) {
        setTxtProgressBar(pb, i)
      }
    }
    if (verbosity > 1) {
      close(pb)
    }
  }
  else {
    if (verbosity > 0) {
      message("Skip calculation of full residual matrix")
    }
    res <- matrix(data = NA, nrow = 0, ncol = 0)
  }
  rv <- list(y = res, model_str = model_str, model_pars = model_pars, 
             model_pars_outliers = model_pars_outliers, model_pars_fit = model_pars_fit, 
             model_str_nonreg = model_str_nonreg, model_pars_nonreg = model_pars_nonreg, 
             arguments = arguments, genes_log_gmean_step1 = genes_log_gmean_step1, 
             cells_step1 = cells_step1, cell_attr = cell_attr)
  rm(res)
  gc(verbose = FALSE)
  times$correct_umi = Sys.time()
  if (return_corrected_umi) {
    if (residual_type != "pearson") {
      message("Will not return corrected UMI because residual type is not set to 'pearson'")
    }
    else {
      rv$umi_corrected <- sctransform::correct(rv, do_round = TRUE, 
                                               do_pos = TRUE, verbosity = verbosity)
      rv$umi_corrected <- as(object = rv$umi_corrected, 
                             Class = "dgCMatrix")
    }
  }
  rv$y[rv$y < res_clip_range[1]] <- res_clip_range[1]
  rv$y[rv$y > res_clip_range[2]] <- res_clip_range[2]
  if (!return_cell_attr) {
    rv[["cell_attr"]] <- NULL
  }
  times$get_gene_attr = Sys.time()
  if (return_gene_attr) {
    if (verbosity > 0) {
      message("Calculating gene attributes")
    }
    gene_attr <- data.frame(detection_rate = genes_cell_count[genes]/ncol(umi), 
                            gmean = 10^genes_log_gmean, variance = sctransform:::row_var(umi))
    if (ncol(rv$y) > 0) {
      gene_attr$residual_mean = rowMeans(rv$y)
      gene_attr$residual_variance = sctransform:::row_var(rv$y)
    }
    if (startsWith(x = method, prefix = "offset")) {
      gene_attr$amean <- rowMeans(umi)
    }
    rv[["gene_attr"]] <- gene_attr
  }
  if (verbosity > 0) {
    message("Wall clock passed: ", capture.output(print(Sys.time() - 
                                                          times$start_time)))
  }
  times$done = Sys.time()
  rv$times <- times
  return(rv)
}

get_model_pars_nonreg_rescale <- function (genes, bin_size, model_pars_fit, regressor_data, umi, 
          model_str_nonreg, cell_attr, verbosity, latent_var_nonreg, min_count = 20) 
{
  bin_ind <- ceiling(x = 1:length(x = genes)/bin_size)
  max_bin <- max(bin_ind)
  if (verbosity > 1) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  model_pars_nonreg <- list()
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- tcrossprod(model_pars_fit[genes_bin, -1, drop = FALSE], 
                     regressor_data)
    umi_bin <- as.matrix(umi[genes_bin, ])
    model_pars_nonreg[[i]] <- do.call(rbind, future.apply:::future_lapply(X = genes_bin, 
                                                           FUN = function(gene) {
                                                             fam <- MASS:::negative.binomial(theta = model_pars_fit[gene, 
                                                                                                             "theta"], link = "log")
                                                             y <- umi_bin[gene, ]
                                                             offs <- mu[gene, ]
                                                             group_ids <- apply(cell_attr[, latent_var_nonreg], 1, which.max)
                                                             group_counts <- aggregate(y, list(group_ids), sum)[,2]
                                                             if(min(group_counts) >= min_count){
                                                               fit <- glm(as.formula(model_str_nonreg), data = cell_attr, 
                                                                          family = fam, offset = offs)
                                                               coeffs <- fit$coefficients
                                                             } else {
                                                               coeffs <- rep(0, length(latent_var_nonreg))
                                                             }
                                                             return(coeffs)
                                                           }, future.seed = TRUE))
    if (verbosity > 1) {
      setTxtProgressBar(pb, i)
    }
  }
  if (verbosity > 1) {
    close(pb)
  }
  model_pars_nonreg <- do.call(rbind, model_pars_nonreg)
  rownames(model_pars_nonreg) <- genes
  return(model_pars_nonreg)
}
