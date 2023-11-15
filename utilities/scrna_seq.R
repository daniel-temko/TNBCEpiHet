#' Summaries average depth and cell numbers across scRNA-seq samples
plot_sc_summaries <- function(sc_list, sample_table, prefix){
  data <- data.frame(sample = factor(sample_table$sample2, levels = sample_table$sample2),
                     batch = sample_table$batch2,
                     mean_depth = sapply(sc_list, function(x) mean(x$nCount_RNA)),
                     n_cells = sapply(sc_list, ncol))
  
  pdf(paste0(prefix, "_depth.pdf"))
  p <- ggplot(data, aes(x = sample, y = mean_depth, fill = batch)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_d3()
  print(p)
  dev.off()
  
  pdf(paste0(prefix, "_cell_num.pdf"))
  p <- ggplot(data, aes(x = sample, y = n_cells, fill = batch)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_d3()
  print(p)
  dev.off()
}
