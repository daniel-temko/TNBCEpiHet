require(foreach)
library(doMC)
registerDoMC(22)

#' Calculate factor scores
call_mofa_factor_scores <- function(call_data){
  weights_mat <- do.call(rbind, call_data$weights)
  exp_mat <- do.call(rbind, call_data$data)
  scores <- factor_scores(weights_mat, exp_mat)
  return(scores)
}

#' Calculates sample x factor matrix from feature x factor weight matrix
#' and feature x sample expression matrix using multiple linear regression.
#' Expression data is feature-centered before the regression -- i.e. 
#' regression is done with no intercept
factor_scores <- function(weights_mat, exp_mat){
  stopifnot(all(rownames(weights_mat) == rownames(exp_mat)))
  lr_mat <- solve(t(weights_mat) %*% weights_mat) %*% t(weights_mat)
  res <- t(lr_mat %*% exp_mat)
  return(res)
}

#' Given a list of data-sets, convert to matrices and remove
#' any sample - column - absent from any dataset, return the
#' list of matrices with columns in the same order across, 
#' data-sets sorted alphabetically
format_model_data_complete_samples <- function(data){
  data <- lapply(data, as.matrix)
  sample_counts <- table(unlist(lapply(data, colnames)))
  complete_samples <- names(sample_counts[which(sample_counts == length(data))])
  data <- lapply(data, function(x) x[,sort(complete_samples)])
  return(data)
}

prepare_mofa_calls <- function(data_list, weight_list){
  stopifnot(all(names(data_list) == names(weight_list)))
  for(i in 1:length(data_list)){
    data_list[[i]] <- sweep(data_list[[i]], 1, rowMeans(data_list[[i]]))
    common_features <- intersect(rownames(data_list[[i]]), rownames(weight_list[[i]]))
    data_list[[i]] <- data_list[[i]][common_features, ]
    weight_list[[i]] <- weight_list[[i]][common_features, ]
  }
  return(list(data = data_list, weights = weight_list))
}

#' Wrapper function to get variance explained for each MOFA factor 
#' in each dataset - following Argelaguet et al. 2018
get_variance_explained_mofa_calls <- function(factor_mat, call_data){
  stopifnot(length(unique(lapply(call_data$weights, colnames))) == 1)
  ve_by_factor <- as.data.frame(lapply(1:length(call_data$data), function(x){
    sapply(1:ncol(call_data$weights[[x]]), function(y){
      var_exp(factor_mat[, y, drop = FALSE], call_data$weights[[x]][, y, drop = FALSE], call_data$data[[x]])
    })  
  }))
  colnames(ve_by_factor) <- names(call_data$data)
  rownames(ve_by_factor) <- colnames(call_data$weights[[1]])
  ve <- sapply(1:length(call_data$data), function(x){
    var_exp(factor_mat, call_data$weights[[x]], call_data$data[[x]])
  })
  names(ve) <- names(call_data$data)
  return(list(r2_total = ve, r2_per_factor = ve_by_factor))
}

#' Calculate variance explained by MOFA factors in a 
#' dataset following Argelaguet 2018. Only consider 
#' features that are present in the expression data
#' for the calculation
var_exp <- function(factor_mat, weights_mat, exp_mat){
  est_vals <- t(factor_mat %*% t(weights_mat))
  rss <- mean((exp_mat - est_vals)^2, na.rm = TRUE)
  tss <- mean(exp_mat^2, na.rm = TRUE)
  return(1 - (rss/tss))
}
