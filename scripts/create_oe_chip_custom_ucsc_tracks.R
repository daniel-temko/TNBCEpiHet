setwd("../")

################################################################################################
# Load data
################################################################################################

load("R_Data/ChIP_PRRX1_OE_Metadata.RData")
load("R_Data/Heatmap_Metadata.RData")

chip_oe_sample_table <- subset(chip_oe_sample_table, sample.type == "ChIP")
chip_oe_sample_table <- chip_oe_sample_table[match(unique(chip_oe_sample_table$experiment.id), chip_oe_sample_table$experiment.id),]
chip_oe_sample_table$tnbc.subtype <- line_annotation[chip_oe_sample_table$parental.cell.line, "type"]

chip_dir <- "../chip_seq_prrx1_oe_cobra"
out_dir <- "analysis/chip_seq_prrx1_oe/custom_tracks"
out_fn <- "all_prrx1_oe_h3k27ac_peaks.bed"

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

stopifnot(!file.exists(paste0(out_dir, "/", out_fn)))
for(i in 1:nrow(chip_oe_sample_table)){
  message(i)
  fn <- with(chip_oe_sample_table, 
             paste0(tolower(sample.id[i]), "/", prefix[i], "_ds0_peaks.narrowPeak"))
  ii <- chip_oe_sample_table$sample.type.code[i]
  sn <- chip_oe_sample_table$experiment.id[i]
  st <- chip_oe_sample_table$tnbc.subtype[i]
  co <- line_colours$type[st]
  rgb <- apply(col2rgb(co), 2, function(x) paste(x, collapse = ","))
  
  pp <- read.table(paste0(chip_dir, "/", fn))
  bb <- genom_df_to_bed(pp)
  
  st <- paste0("track name=", ii, " ", 
               "description=", paste0('"', sn, " peaks", '"'), " ", 
               "visibility=", 1, " ",
               "color=", rgb)
  
  write.table(st, paste0(out_dir, "/", out_fn), 
              quote = FALSE, 
              col.names = FALSE, 
              row.names = FALSE,
              append = TRUE)
  
  write_bed(bb, paste0(out_dir, "/", out_fn), append = TRUE)
}