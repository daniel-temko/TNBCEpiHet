#' Aggregate binary peak matrix to the closest gene column. If canon
#' is present filter before aggregation. Ignore peaks with gene.symbol
#' "."
agg_peaks_to_genes <- function(peak_df, ann_df, canon = NULL){
  if(!is.null(canon)){
    keep <- ann_df$chr.region %in% canon
    peak_df <- peak_df[keep, ]
    ann_df <- ann_df[keep, ]
  }
  peak_df$closest.gene <- ann_df$gene.symbol
  df <- aggregate(. ~ closest.gene, data = peak_df, FUN = function(x) 1 - prod(1 - x))
  if("." %in% df$closest.gene){
    df <- df[-match(".", df$closest.gene), ]  
  }
  rownames(df) <- df[,1]
  df[,1] <- NULL
  return(as.data.frame(df > 0))
}

#' Add an an annotation column to targets_df containing the
#' cell lines for which each region is a target
annotate_targets <- function(targets_df, cell_lines){
  print("Annotating rows with cell lines...")
  line_targets <- lapply(1:nrow(targets_df), function(x){
    sort(unique(cell_lines[which(targets_df[x,] == 1)]))
  })
  num_lines <- sapply(line_targets, length)
  line_str <- sapply(line_targets, function(x) paste(x, collapse = ":"))
  targets_df$targets <- line_str
  targets_df$num_lines <- num_lines
  return(targets_df)
}

#' Add expression information to targets_df ChIP-seq results
add_expression <- function(targets_df, exp_df, col_name){
  targets_df[[col_name]] <- "ns"
  exp_df$dir <- sapply(exp_df$log2FoldChange, function(x){
    if(x > 0) {
      "up"
    } else {
      "down"
    }
  })
  common_ids <- intersect(rownames(targets_df), rownames(exp_df))
  targets_df[common_ids, col_name] <- exp_df[common_ids, "dir"]
  return(targets_df)
}