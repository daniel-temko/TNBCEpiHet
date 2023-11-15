#' Extract FDR val for each sig in sig_ids, output
#' results as ordered numeric vector. Output NA for
#' any signature absent from the df, and "No tested"
#' for signatures with too few or too many genes
get_fdrs_from_df <- function(in_df, sig_ids, rejected_sig_ids){
  res <- sapply(sig_ids, function(x){
    if(x %in% rejected_sig_ids){
      "Not tested"
    } else if (!(x %in% in_df$NAME)){
      # Not enriched
      NA
    } else {
      in_df$FDR.q.val[match(x, in_df$NAME)]
    }
  })
  names(res) <- sig_ids
  return(res)
}

#' Identify the signatures that were not tested from
#' gene_set_df
get_rejected_sig_ids <- function(gene_set_df, sig_ids){
  stopifnot(all(sort(gene_set_df$NAME) == sort(sig_ids)))
  statuses <- gene_set_df$STATUS[match(sig_ids, gene_set_df$NAME)]
  rejected <- which(statuses == "Rejected!")
  return(sig_ids[rejected])
}

#' For each row in contrasts df, create two corresponding columns 
#' - one per comparison direction - in an output df containing the FDR values 
#' for each of the investigated signatures in the rows
get_fdrs <- function(contrast_df, sig_names, gsea_path, direction_id){
  sig_ids <- toupper(sig_names)
  res <- lapply(1:nrow(contrast_df), function(x){
    dir_id <- paste0(contrast_df$contrast_id[x], ".Gsea.", contrast_df$time.stamp[x])
    gene_set_df <- read.table(file.path(gsea_path, dir_id, "gene_set_sizes.tsv"), sep = "\t", header = TRUE)
    rejected_sig_ids <- get_rejected_sig_ids(gene_set_df, sig_ids)
    fdr_fname <- paste0("gsea_report_for_", contrast_df[[direction_id]][x], "_", contrast_df$time.stamp[x], ".tsv")
    fdr_df <- read.table(file.path(gsea_path, dir_id, fdr_fname), sep = "\t", header = TRUE)
    get_fdrs_from_df(fdr_df, sig_ids, rejected_sig_ids)
  })
  df <- as.data.frame(res)
  colnames(df) <- paste0(contrast_df$cell.line, "_", contrast_df[[paste0("treatment_", direction_id)]][1])
  stopifnot(all(rownames(df) == sig_ids))
  rownames(df) <- sig_names
  return(df)
}

#' Convert a data frame containing FDR values to one
#' containing TRUE/FALSE indicators of significance
get_significant <- function(fdrs, sig_val = 0.05){
  res <- fdrs < sig_val
  nt_ids <- which(fdrs == "Not tested")
  res[nt_ids] <- "Not tested"
  return(res)
}

#' Create GSEA plot from data frame
plot_gsea_from_df <- function(contrast_df, sig_df, sig_id, fdr_direction, fdr_val, line_col = "blue"){
  x_min <- min(sig_df$RANK.IN.GENE.LIST)
  x_max <- max(sig_df$RANK.IN.GENE.LIST)
  y_min <- min(sig_df$RUNNING.ES)
  y_max <- max(sig_df$RUNNING.ES) 
  fdr_str <- paste0("Q=", format(fdr_val, scientific = TRUE, digits = 3))
  plot(sig_df$RANK.IN.GENE.LIST, sig_df$RUNNING.ES, type = "l", ylim = c(y_min, y_max), ylab = "", xlab = "", col = line_col)
  abline(h = 0, col = "grey")
  segments(x0 = sig_df$RANK.IN.GENE.LIST, y0 =  y_min - 1, x1 = sig_df$RANK.IN.GENE.LIST, y1 = y_min, lwd = 0.5)
  if(fdr_direction == "a"){
    text(x = x_min + (x_max - x_min) * 0.8, y = y_max * 0.8, fdr_str)
  } else {
    text(x = x_min + (x_max - x_min) * 0.2, y = y_min * 0.8, fdr_str)
  }
  title(ylab = "Running ES")
  title(main = paste0("Enrichment ", tolower(sig_id), " signature"))
  title(xlab = paste0("Gene ranked from up to down\nregulation (", contrast_df$a, "/", contrast_df$b, ")"))
}

#' Crate GSEA plot 
plot_gsea <- function(contrast_df, sig_id, gsea_path){
  message("Plotting...")
  stopifnot(nrow(contrast_df) == 1)
  dir_id <- paste0(contrast_df$contrast_id, ".Gsea.", contrast_df$time.stamp)
  sig_fname <- paste0(sig_id, ".tsv")
  sig_df <- read.table(file.path(gsea_path, dir_id, sig_fname), sep = "\t", header = TRUE, row.names = 1)
  
  fdr_fname_a <- paste0("gsea_report_for_", contrast_df$a, "_", contrast_df$time.stamp, ".tsv")
  fdr_df_a <- read.table(file.path(gsea_path, dir_id, fdr_fname_a), sep = "\t", header = TRUE)
  fdr_a <- fdr_df_a$FDR.q.val[match(sig_id, fdr_df_a$NAME)]
  
  fdr_fname_b <- paste0("gsea_report_for_", contrast_df$b, "_", contrast_df$time.stamp, ".tsv")
  fdr_df_b <- read.table(file.path(gsea_path, dir_id, fdr_fname_b), sep = "\t", header = TRUE)
  fdr_b <- fdr_df_b$FDR.q.val[match(sig_id, fdr_df_b$NAME)]
  
  fdr_vals <- c(a = fdr_a, b = fdr_b)
  stopifnot(length(which(is.na(fdr_vals))) == 1)
  fdr_val <- fdr_vals[which(!is.na(fdr_vals))]
  fdr_direction <- names(fdr_val)
  
  plot_gsea_from_df(contrast_df, sig_df, sig_id, fdr_direction, fdr_val)
}