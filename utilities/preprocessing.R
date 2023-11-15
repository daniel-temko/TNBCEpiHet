library(reshape2)
library(ggplot2)

#' Select HVG's by variance
variable_features <- function(feat_data){
  feat_var <- apply(feat_data, 1, var)
  names(sort(feat_var, decreasing = TRUE))
}

#' Get the first perc_features features from vector of 
#' features 'features'
first_n_percent <- function(features, perc_features){
  num_features <- round(perc_features * (1/100) * length(features))
  return(features[1:num_features])
}

#' Filter in_data to the top perc_features most variable features
#' and rank the rows in the output by variance
get_heatmap_data <- function(in_data, perc_features = 20){
  ranked_features <- variable_features(in_data)
  top_features <- first_n_percent(ranked_features, perc_features)
  return(in_data[top_features, ])
}

#' Plot total and mapped reads
PlotMappedReads <- function(samples, sample_table, data_dir, analysis_dir){
  chip_total_reads <- c()
  chip_mapped_reads <- c()
  for(i in 1:length(samples)){
    print(i)
    sample_info <- subset(sample_table, Sample.Name == samples[i])
    file_name <- with(sample_info, paste0(tolower(Sample.Name[1]), "/", prefix[1], ".rmdup.mapstat"))
    map_info <- read.table(paste0(data_dir, "/", file_name), fill = TRUE, stringsAsFactors = FALSE)
    chip_total_reads <- c(chip_total_reads, as.numeric(map_info[1,1]))
    chip_mapped_reads <- c(chip_mapped_reads, as.numeric(map_info[5,1]))
  }
  
  chip_mapping_info <- data.frame(sample = samples,
                                  total = chip_total_reads,
                                  mapped = chip_mapped_reads)
  
  Data <- melt(chip_mapping_info, id.vars = "sample")
  Data$sample <- factor(Data$sample, levels = samples)
  
  pdf(file = paste0(analysis_dir, "mapping_rate.pdf"), width = 11, height = 8)
  p <- ggplot(data = Data, aes(x = sample, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme(axis.text.x = element_text(angle = 90))
  print(p)
  dev.off()
}

#' Plot duplicate read rate
PlotDupRate <- function(samples, sample_table, data_dir, analysis_dir){
  dup_rate <- c()
  for(i in 1:length(samples)){
    print(i)
    sample_info <- subset(sample_table, Sample.Name == samples[i])
    file_name <- with(sample_info, paste0(tolower(Sample.Name[1]), "/", prefix[1], ".dup.txt"))
    dup_info <- read.table(paste0(data_dir, "/", file_name), 
                           stringsAsFactors = FALSE,
                           header = TRUE,
                           fill = TRUE,
                           sep = "\t")
    dup_rate <- c(dup_rate, dup_info$PERCENT_DUPLICATION[1])
  }
  
  chip_dup_info <- data.frame(sample = samples, dup.rate = dup_rate)
  Data <- melt(chip_dup_info, id.vars = "sample")
  Data$sample <- factor(Data$sample, levels = samples)
  
  pdf(file = paste0(analysis_dir, "dup_rate.pdf"), width = 11, height = 8)
  p <- ggplot(data = Data, aes(x = sample, y = value)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme(axis.text.x = element_text(angle = 90))
  print(p)
  dev.off()
}
