setwd("../")

################################################################################################
# Load data
################################################################################################

load("R_Data/PRRX1_ChIP_Metadata.RData")
prrx1_chip_sample_table <- subset(prrx1_chip_sample_table, (sample.type == "ChIP") & (sonication.condition %in% c(1, 2)))
prrx1_chip_sample_table$sample.id2 <- with(prrx1_chip_sample_table, paste0(cell.line, "_", sonication.condition))

colour_code <- read.table("metadata/colour_code_coclustering.csv", sep = ",", header = TRUE, comment.char = "")
colour_map <- as.character(colour_code$colour)
names(colour_map) <- colour_code$type
names(colour_map)[1:3] <- c("basal", "mesenchymal", "luminal")

chip_dir <- "../chip_seq_prrx1_cobra"
out_dir <- "analysis/prrx1_chip_seq/custom_tracks"
out_fn <- "all_prrx1_peaks.bed"

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

oo <- order(prrx1_chip_sample_table$cell.line.type, prrx1_chip_sample_table$cell.line)
prrx1_chip_sample_table <- prrx1_chip_sample_table[oo,]

stopifnot(!file.exists(paste0(out_dir, "/", out_fn)))
for(i in 1:nrow(prrx1_chip_sample_table)){
  message(i)
  fn <- with(prrx1_chip_sample_table, 
             paste0(sample.id[i], "/", prefix[i], "_ds0_peaks.narrowPeak"))
  sn <- prrx1_chip_sample_table$sample.id2[i]
  tt <- prrx1_chip_sample_table$cell.line.type[i]
  co <- colour_map[tt]
  rgb <- apply(col2rgb(co), 2, function(x) paste(x, collapse = ","))
  
  pp <- read.table(paste0(chip_dir, "/", fn))
  bb <- genom_df_to_bed(pp)
  
  st <- paste0("track name=", sn, " ", 
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
