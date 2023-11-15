library(doMC)
library(GenomicRanges)
library(varhandle)
library(ggplot2)
registerDoMC(22)

#' Create path if it does not exist
create_path <- function(path_name, recursive = TRUE){
  if(!dir.exists(path_name)) {
    message("Creating ", path_name)
    dir.create(path_name, recursive = recursive)
  }
  return(path_name)
}

#' Substitute into character vectors - throws an error if there is an
#' unhandled case
easy_find_replace <- function(in_vec, find_vec, replace_vec, na_val = NA){
  stopifnot(is.character(in_vec))
  stopifnot(is.character(find_vec))
  stopifnot(length(unique(find_vec)) == length(find_vec))
  sapply(in_vec, function(x) {
    if(is.na(x)){
      na_val
    } else if(x %in% find_vec){
      replace_vec[match(x, find_vec)]
    } else {
      stop(paste0("Error - unhandled case: ", x))
    }
  })
}

#' Read CSV file
read_csv <- function(file){
  read.table(file,
             sep = ",",
             header = TRUE,
             stringsAsFactors = FALSE)
}

#' Read MSIGDB signatures
read_msigdb <- function(fn){
  n_col <- max(count.fields(fn, sep = "\t"))
  read.table(fn, sep = "\t", fill = TRUE, col.names = paste0("C", 1:n_col))
}

#' Output file to csv
write_csv <- function(df, file_name){
  write.table(df, file_name, sep = ",", row.names = FALSE, quote=FALSE,)
}

#' Parse MSIGDB-formatted gene sets
parse_msigdb <- function(gs_mat, all_genes){
  stopifnot(max(table(gs_mat[,1])) == 1)
  stopifnot(class(gs_mat) == c("matrix", "array"))
  res <- foreach(i = 1:nrow(gs_mat)) %dopar% {
    genes <- as.character(gs_mat[i, -(1:2)])
    genes <- genes[which(genes != "")]
    set_name <- as.character(gs_mat[i,1])
    data.frame(set = set_name, genes = genes)
  }
  res_df <- do.call(rbind, res)
  ann_genes <- sort(unique(c(res_df$genes, all_genes)))
  res_df$genes <- factor(res_df$genes, levels = ann_genes)
  res_df_mat <- as.data.frame.matrix(table(res_df$set, res_df$genes))
  as.matrix(res_df_mat)
}

#' Bring column to front of data.frame
bring_to_front <- function(df, col_name){
  df <- df[,c(col_name, colnames(df)[which(colnames(df) != col_name)])]
  return(df)
}

#' Bring mutliple columns to front
bring_multiple_to_front <- function(df, col_names){
  df <- df[,c(col_names, colnames(df)[which(!(colnames(df) %in% col_names))])]
  return(df)
}

#' Remove duplicate data frame rows
dedup_df <- function(df, col_num=1){
  df <- df[match(unique(df[, col_num]), df[, col_num]), , drop=FALSE]
}

#' Rename data frame column
rename_col <- function(df, col_num, col_name){
  colnames(df)[col_num] <- col_name
  return(df)
}

#' Assign row names from data frame column
rownames_from_col <- function(df, col_num = 1){
  rownames(df) <- df[,col_num]
  return(df)
}

#' Concatenate 'chr', 'start', and 'end' into a region id, without 
#' invoking scientific formatting
region_ids_from_zero_based_half_open <- function(chr, start, end){
  if(!is.character(chr)){
    chr <- as.character(chr)
  }
  if(!is.integer(start)){
    start <- as.integer(start)
  }
  if(!is.integer(end)){
    end <- as.integer(end)
  }
  return(sprintf("%s:%d-%d", chr, start, end))
}

#' Subtract one from start, to convert to zero-based half open.
#' Concatenate 'chr', 'start', and 'end' into a region id, without 
#' invoking scientific formatting
region_ids <- function(chr, start, end){
  start <- start - 1
  if(!is.character(chr)){
    chr <- as.character(chr)
  }
  if(!is.integer(start)){
    start <- as.integer(start)
  }
  if(!is.integer(end)){
    end <- as.integer(end)
  }
  return(sprintf("%s:%d-%d", chr, start, end))
}

#' Create to UCSF format 1-based region ids
pub_region_ids <- function(chr, start, end){
  if(!is.character(chr)){
    chr <- as.character(chr)
  }
  if(!is.integer(start)){
    start <- as.integer(start)
  }
  if(!is.integer(end)){
    end <- as.integer(end)
  }
  return(sprintf("%s:%d-%d", chr, start, end))
}

#' Get region ID's from GRanges, conversion to 0-based start,
#' is done within the embedded call to region_ids
region_ids_from_granges <- function(gr){
  chr <- as.character(seqnames(gr))
  start <- start(gr)
  end <- end(gr)
  return(region_ids(chr, start, end))
}

#' Split region ids into chromosome, start, and end vectors
#' Use these components to construct GRanges, with starts incremented
#' by one to convert to 1-based
granges_from_region_ids <- function(in_ids){
  components <- lapply(in_ids, function(x) strsplit(x, ":")[[1]])
  chr <- sapply(components, function(x) x[1])
  pos_components <- lapply(components, function(x) strsplit(x[2], "-")[[1]])
  start <- as.integer(sapply(pos_components, function(x) x[1]))
  end <- as.integer(sapply(pos_components, function(x) x[2]))
  in_df <- data.frame(chr = chr, start = start, end = end, row.names = NULL)
  gr <- makeGRangesFromDataFrame(in_df, starts.in.df.are.0based = TRUE)
  return(gr)
}

#' Convert internal region ID's to pub region
#' ID's
pub_region_ids_from_region_ids <- function(in_ids){
  components <- lapply(in_ids, function(x) strsplit(x, ":")[[1]])
  chr <- sapply(components, function(x) x[1])
  pos_components <- lapply(components, function(x) strsplit(x[2], "-")[[1]])
  start <- as.integer(sapply(pos_components, function(x) x[1]))
  end <- as.integer(sapply(pos_components, function(x) x[2]))
  start <- start + 1
  return(pub_region_ids(chr = chr, start = start, end = end))
}

#' Return quantile cuts for values
get_quantiles <- function(values, n_quantiles){
  breaks <- quantile(values, probs = seq(0, 1, 1/n_quantiles))
  cuts <- cut(values, breaks = breaks, include.lowest = TRUE)
  if(n_quantiles > 2){
    out_vals <- factor(cuts, labels = paste0("q", 1:n_quantiles))  
  } else{
    out_vals <- factor(cuts, labels = c("low", "high"))
  }
  return(out_vals)
}

#' Removes rows with any missing values
#' and drops unused levels
drop_missing <- function(in_data){
  stopifnot(is.data.frame(in_data))
  out_data <- droplevels(in_data[complete.cases(in_data), ])
  return(out_data)
}

#' Confirm that columns status_id and time_id are numeric. Confirm at least one entry
#' in var_ids and that all columns referenced in var_ids and cov_ids are factor or 
#' numeric type
validate_survival_data <- function(status_id, time_id, var_ids, in_data, cov_ids = c()){
  stopifnot(all(c(status_id, time_id, var_ids, cov_ids) %in% colnames(in_data)))
  stopifnot(is.numeric(in_data[, status_id]))
  stopifnot(is.numeric(in_data[, time_id]))
  stopifnot(length(var_ids) >= 1)
  stopifnot(all(sapply(c(var_ids, cov_ids), function(x) is.numeric(in_data[, x]) | is.factor(in_data[, x]))))
}

#' Summarize results of CoxPH model using only numeric and factor variaables - probably will
#' not work with splines
coxph_summary <- function(cox_res){
  coef_cols <- c("exp(coef)", "lower .95", "upper .95")
  res <- cbind(summary(cox_res)$conf.int[, coef_cols, drop = FALSE], 
               summary(cox_res)$coefficients[, "Pr(>|z|)", drop = FALSE])
  res <- as.data.frame(res)
  return(res)
}

#' Rank values in vec, in descending order
get_rank <- function(vec){
  length(vec) + 1 - rank(vec)
}

#' Calculate vector of all pairwise Euclidean distances involving different 
#' cells
pairwise_distances_parallel <- function(df, n_threads = 20){
  require(foreach)
  require(doMC)
  registerDoMC(n_threads)
  res = foreach(i = 1:(ncol(df) - 1)) %dopar% {
    diff <- df[, (i + 1):ncol(df), drop = FALSE] - df[, i]
    pw_distances <- apply(diff, 2, function(x) sqrt(sum(x^2)))
  }
  unlist(res)
}

#' Optimized version of pairwise distances function
pairwise_distances_optim <- function(df, verbose = FALSE){
  res <- c()
  for(i in 1:(ncol(df) - 1)){
    if(verbose) {
      print(paste0(i, " of ", (ncol(df) -1)))
    } else {
      if((i %% 100) == 0) print(paste0(i, " of ", ncol(df) - 1))
    }
    diff <- df[, (i+1):ncol(df), drop = FALSE] - df[, i]
    pw_distances <- sqrt(colSums(diff^2))
    res <- c(res, pw_distances)
  }
  return(res)
}

#' Compute vector of all pairwise distances
pairwise_distances <- function(df, verbose = FALSE){
  res <- c()
  for(i in 1:(ncol(df) - 1)){
    if(verbose) {
      print(paste0(i, " of ", (ncol(df) -1)))
    } else {
      if((i %% 100) == 0) print(paste0(i, " of ", ncol(df) - 1))
    }
    diff <- df[, (i+1):ncol(df), drop = FALSE] - df[, i]
    pw_distances <- apply(diff, 2, function(x) sqrt(sum(x^2)))
    res <- c(res, pw_distances)
  }
  return(res)
}

#' Create sbatch job submission header with given specifications
#' 
#' Note: mem is specificied in GB!
#' Note: requires directory 'logs' to exist
slurm_header <- function(job_name, n_tasks, mem, days = 1, mem_type = "mpc", exclude=NULL){
  if(mem_type == "mpc"){
    mem_str <- "#SBATCH --mem-per-cpu="
  } else if (mem_type == "mem") {
    mem_str <- "#SBATCH --mem="
  } else {
    warning("Unrecognised memtype...")
  }
  header <- paste0("#!/bin/bash\n",
                   "#\n",
                   "#SBATCH --ntasks=", n_tasks, "\n",
                   "#SBATCH --time=", days, "-00:00\n",
                   mem_str, mem, "000\n",
                   "#SBATCH --job-name=", job_name, "\n",
                   "#SBATCH -e logs/", job_name, "_%j.err\n",
                   "#SBATCH -o logs/", job_name, "_%j.out\n")
  if(!is.null(exclude)){
    header <- paste0(header,
                     "#SBATCH -x ", paste(exclude, collapse=","), "\n")
  }
  header <- paste0(header, "\n")
  return(header)
}

#' Get first n elements of vector
get_topn <- function(in_vec, n){
  if(n == 0){
    return(in_vec[c()])
  } else {
    return(in_vec[1:n])
  }
}

#' Restricts gsl to elements in the genes vector. If
#' n_targets is given, then takes the first n_target entries
#' that are in the genes vectors, or all such entries if n_target >
#' length(genes)
subset_gsl <- function(gsl, genes, n_target = NA){
  for(i in 1:length(gsl)){
    for(j in 1:length(gsl[[i]])){
      vec <- gsl[[i]][[j]]
      new_vec <- vec[which(vec %in% genes)]
      if(!is.na(n_target)){
        new_vec <- get_topn(new_vec, min(n_target, length(new_vec)))
      }
      gsl[[i]][[j]] <- new_vec
    }
  }
  return(gsl)
}

#' Calculates the minimum length of overlap between each gene vector in each gene set 'n'
#' and selects the top n such elements from each gene vector in each gene set
uniform_subset_gsl <- function(gsl, genes){
  cts <- sapply(gsl, function(x){
    sapply(x, function(y){
      length(intersect(genes, y))
    })
  })
  sel_num <- min(cts)
  subset_gsl(gsl, genes, n_target = sel_num)
}

#' Append (N) onto names of duplicates features to make
#' them unique
rename_duplicates <- function(vec){
  dups <- table(vec)[which(table(vec) > 1)]
  if(length(dups) == 0) return(vec)
  for(i in 1:length(dups)){
    ids <- which(vec == names(dups)[i])
    vec[ids] <- paste0(vec[ids], " (", 1:dups[i], ")")
  }
  return(vec)
}

#' regress out the effect of covariate from each row
regress_matrix_covariate <- function(in_data, covariate){
  mean_vals <- rowMeans(in_data) 
  resid_mat <- t(apply(in_data, 1, function(x) lm(x ~ covariate)$residuals))
  resid_mat + mean_vals
}

#' Checks that all columns in in_data match by name to an 
#' element in covariate - throws an error if not. Reorders
#' elements by name and regresses out the effect of covariate
#' from each row
#' 
#' Note: Requires matrix/data.frame has column names and 
#' covariates is a named vector
regress_named_matrix_covariate <- function(in_data, covariate){
  stopifnot(all(colnames(in_data) %in% names(covariate)))
  reordered_covariate <- covariate[colnames(in_data)]
  regress_matrix_covariate(in_data, reordered_covariate)
}

#' Find order of sample groups vector
group_order <- function(groups, ordered_samples){
  nsamples_padded <- sapply(groups, 
                            function(x) str_pad(str_count(as.character(x), ":") + 1, 2, side = "left", pad = "0"))
  combn_id <- sapply(groups, 
                     function(x) paste(str_pad(match(strsplit(as.character(x), ":")[[1]], ordered_samples), 2, side = "left", pad = "0"), collapse = "_"))
  return(order(paste0(nsamples_padded, combn_id)))
}

#' Make column names compatible with R processing
r_process_col_names <- function(col_names){
  new_names <- sapply(col_names, function(x){
    if(check.numeric(substr(x, 1, 1))){
      paste0("X", x)
    } else
      x
  })
  sapply(new_names, function(x) gsub("(\\(|\\)|-| |/)", ".", x))
}

#' Get ids of repeated elements in vec
get_tie_ids <- function(vec){
  tie_vals <- names(table(vec))[which(table(vec) > 1)]
  which(vec %in% tie_vals)
}

#' Retrieve the indexes from vector with the n_top and n_bottom
#' top and bottom values respectively
top_and_bottom_ids <- function(vec, n_top = 10, n_bottom = 10){
  stopifnot((n_top + n_bottom) <= length(vec))
  ord <- order(vec, decreasing = TRUE)
  if(n_top > 0){
    top_ids <- ord[1:n_top]  
  } else {
    top_ids <- c()
  }
  if(n_bottom > 0){
    bottom_ids <- rev(ord)[1:n_bottom]  
  } else {
    bottom_ids <- c()
  }
  return(list(top_ids = top_ids, bottom_ids = bottom_ids))
}

#' ggplot2 theme
my_theme <- function(text_size = 20){
  theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(size = text_size))
}

my_theme2 <- function (base_size = 11, base_family = "", base_line_size = base_size/22, 
                      base_rect_size = base_size/22) {
  half_line <- base_size/2
  theme_bw(base_size = base_size, base_family = base_family, 
           base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(axis.text = element_text(colour = "black", size = rel(0.8)), 
          axis.ticks = element_line(colour = "black", size = rel(0.5)), 
          panel.border = element_rect(fill = NA, colour = "black", 
                                      size = rel(1)), panel.grid = element_blank(),
          strip.background = element_rect(fill = "black"), 
          strip.text = element_text(colour = "white", size = rel(0.8), 
                                    margin = margin(0.8 * half_line, 0.8 * half_line, 
                                                    0.8 * half_line, 0.8 * half_line)), complete = TRUE)
}

#' Convert png file to pdf
png2pdf <- function(png_fname, pdf_fname){
  img <- png::readPNG(png_fname)
  pdf(pdf_fname)
  print(ggplot(data.frame()) + theme_void() + annotation_raster(img, -Inf, Inf, -Inf, Inf))
  dev.off()
}

#' Apply rule of thumb from for when Chi-squared test can be used on
#' contingency table
validate_chisq <- function(cont_mat){
  rm_props <- rowSums(cont_mat) / sum(cont_mat)
  cm_props <- colSums(cont_mat) / sum(cont_mat)
  exp_mat <- matrix(sum(cont_mat), nrow = nrow(cont_mat), ncol = ncol(cont_mat))
  exp_mat <- exp_mat * rm_props
  exp_mat <- t(t(exp_mat) * cm_props)
  res <- (length(which(exp_mat < 5)) / length(exp_mat)) <= 0.2
  list(exp_mat, res)
}

#' Randomly down-sample a vector to a given length, leaving
#' the order unchanged
down_sample <- function(vec, target_len, rseed = 1234){
  set.seed(rseed)
  ii <- sample(1:length(vec), min(length(vec), target_len))
  return(vec[sort(ii)])
}

#' Intersect vec1 and vec2 and return the intersection
#' in order of average position of each element in the two
#' vectors. Ties in average position are broken randomly
ranked_intersect <- function(vec1, vec2, rseed = 1234){
  set.seed(rseed)
  ii <- intersect(vec1, vec2)
  s1 <- vec1[which(vec1 %in% ii)]
  s2 <- vec2[which(vec2 %in% ii)]
  r1 <- match(ii, s1)
  r2 <- match(ii, s2)
  aa <- (r1 + r2) / 2
  return(ii[order(rank(aa, ties.method = "random"))])
}

#' Find union of vec1 and vec2 and return the result
#' sorted by the average position of each element in
#' the two vectors, where applicable, or the position
#' of each element in the vector where it occurs
#' otherwise. Ties are broken randomly
ranked_union <- function(vec1, vec2, rseed = 1234){
  set.seed(rseed)
  uu <- union(vec1, vec2)
  r1 <- match(uu, vec1)
  r2 <- match(uu, vec2)
  aa <- rowMeans(data.frame(r1, r2), na.rm = TRUE)
  return(uu[order(rank(aa, ties.method = "random"))])
}

