#' Convert to 0-based start, and output to tab-delimited file
#' suppressing scientific notation. 
#' Note: 
#' 1. There is a risk of integer overflow if dealing with large 
#'    non-human genomes. In this case throws an overflow error
#' 2. Relies on scientific formatting suppression in read.table
#'    for integer type.
write_df_to_bed <- function(in_df, file_path, starts_in_df_are_1based = TRUE){
  if(starts_in_df_are_1based){
    in_df[,2] <- in_df[,2] - 1
  }
  if (!is.integer(in_df[,2])) {
    in_df[,2] <- as.integer(in_df[,2])
  }
  if(!is.integer(in_df[,3])){
    in_df[,3] <- as.integer(in_df[,3])  
  }
  write.table(in_df,
              file = file_path,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t")
}