setwd("../")

################################################################################################
# Load data
################################################################################################

load("R_Data/ChIP_Metadata.RData")
load("R_Data/Heatmap_Metadata.RData")
chip_sample_table <- subset(chip_sample_table, Pulldown == "H3K27ac")

chip_dir <- "../chip_seq_cell_lines_cobra"
out_dir <- "analysis/chip_seq/custom_tracks"
out_fn <- "all_h3k27ac_peaks.bed"

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

oo <- order(chip_sample_table$tnbc.subtype, chip_sample_table$Cell.Line)
chip_sample_table <- chip_sample_table[oo,]

stopifnot(!file.exists(paste0(out_dir, "/", out_fn)))
for(i in 1:nrow(chip_sample_table)){
  message(i)
  fn <- with(chip_sample_table, 
             paste0(tolower(Sample.Name[i]), "/", prefix[i], "_ds0_peaks.narrowPeak"))
  cl <- chip_sample_table$Cell.Line[i]
  st <- chip_sample_table$tnbc.subtype[i]
  co <- line_colours$type[st]
  rgb <- apply(col2rgb(co), 2, function(x) paste(x, collapse = ","))
  
  pp <- read.table(paste0(chip_dir, "/", fn))
  bb <- genom_df_to_bed(pp)
  
  st <- paste0("track name=", cl, " ", 
               "description=", paste0('"', cl, " peaks", '"'), " ", 
               "visibility=", 1, " ",
               "color=", rgb)
  
  write.table(st, paste0(out_dir, "/", out_fn), 
              quote = FALSE, 
              col.names = FALSE, 
              row.names = FALSE,
              append = TRUE)
  
  write_bed(bb, paste0(out_dir, "/", out_fn), append = TRUE)
}
