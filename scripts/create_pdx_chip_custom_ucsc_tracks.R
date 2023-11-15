setwd("../")

################################################################################################
# Load data
################################################################################################

load("R_Data/PDX_ChIP_Metadata.RData")
pdx_chip_sample_table <- subset(pdx_chip_sample_table, Pulldown == "H3K27ac")
pdx_chip_sample_table$sample.id <- sapply(pdx_chip_sample_table$PDX.ID, function(x) gsub("-", "", x))
rgb <- "130,127,127"

chip_dir <- "../pdx_chip_seq_cobra"
out_dir <- "analysis/pdx_chip_seq/custom_tracks"
out_fn <- "all_pdx_h3k27ac_peaks.bed"

################################################################################################
# Local functions
################################################################################################

genom_df_to_bed <- function(df){
  chr_vals <- as.character(df[,1])
  start_pos <- as.integer(df[,2])
  end_pos <- as.integer(df[,3])
  id <- as.character(df[,4])
  df <- data.frame(chr = chr_vals, start = start_pos, end = end_pos, id = id)
  return(df)
}

write_bed <- function(df, file_name, append=FALSE){
  write.table(df,
              file_name, 
              sep = "\t", 
              quote = FALSE,
              col.names = FALSE,
              row.names = FALSE,
              append = append)
}

################################################################################################
# Transform data
################################################################################################

oo <- order(pdx_chip_sample_table$PDX.ID)
pdx_chip_sample_table <- pdx_chip_sample_table[oo, ]

stopifnot(!file.exists(paste0(out_dir, "/", out_fn)))
for(i in 1:nrow(pdx_chip_sample_table)){
  message(i)
  fn <- with(pdx_chip_sample_table, 
             paste0(Sample.Name[i], "/", prefix[i], "_ds0_peaks.narrowPeak"))
  ii <- pdx_chip_sample_table$sample.id[i]
  
  pp <- read.table(paste0(chip_dir, "/", fn))
  bb <- genom_df_to_bed(pp)
  
  st <- paste0("track name=", ii, " ", 
               "description=", paste0('"', ii, " peaks", '"'), " ", 
               "visibility=", 1, " ",
               "color=", rgb)
  
  write.table(st, paste0(out_dir, "/", out_fn), 
              quote = FALSE, 
              col.names = FALSE, 
              row.names = FALSE,
              append = TRUE)
  
  write_bed(bb, paste0(out_dir, "/", out_fn), append = TRUE)
}


