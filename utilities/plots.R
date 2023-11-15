library(ggplot2)

#' plot Mclust results
plot_mclust_results <- function(res){
  plot(res, what = "BIC")
  plot(res, what = "classification")
  plot(res, what = "density")
}

#' Boxplot vlaues by type, with colouring by type
plot_by_group <- function(vals, 
                          types,
                          ylab,
                          cols = line_colours$type[c("luminal", "basal", "mesenchymal")]){
  data <- data.frame(val = vals, type = types)
  y_min <- min(vals)
  y_max <- max(vals)
  data$type <- factor(data$type, levels = names(cols))
  par(mar = c(6, 8, 4, 4))
  boxplot(val ~ type, data = data, col = cols, ylab = ylab, cex.lab = 2, cex.axis = 2, outline = FALSE,
          ylim = c(y_min, y_max))
  stripchart(val ~ type, data = data, method = "jitter", pch = 19, vertical = TRUE, add = TRUE)
}

#' Barplot values by cell line, colouring by type
plot_by_line <- function(vals, 
                         types, 
                         ylab,
                         samp_names,
                         cols = line_colours$type[c("luminal", "basal", "mesenchymal")]){
  data <- data.frame(val = vals, type = types, sample = samp_names)
  data <- data[order(data$val),]
  par(mar = c(6, 4, 4, 4))
  barplot(data$val, names.arg = data$sample, las = 2, col = cols[as.character(data$type)], ylab = ylab)
}

plot_violin_by_group <- function(vals,
                                 types,
                                 ylab,
                                 p_adj,
                                 cols = line_colours$type[c("luminal", "basal", "mesenchymal")]){
  data <- data.frame(val = vals, type = types)
  data$type <- factor(data$type, levels = names(cols))
  y_min <- min(data$val)
  y_max <- max(data$val)
  plot_n_ann <- data.frame(x = 1:length(levels(data$type)),
                           y = y_min + (y_max - y_min) * 1.2,
                           type = levels(data$type),
                           label = paste0("N=", as.numeric(table(data$type))))
  plot_p_ann <- data.frame(x = median(1:length(levels(data$type))),
                           y = y_min + (y_max - y_min) * 1.4,
                           type = NA,
                           label = paste0("P.adj=", format(p_adj, scientific = TRUE, digits = 3)))
  p <- ggplot(data, aes(x = type, y = val, fill = type)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, show.legend = FALSE) +
    scale_fill_manual(values = cols) +
    ylab(ylab) +
    geom_text(size = 6,
              data = rbind(plot_n_ann, plot_p_ann),
              aes(x = x, y = y, label = label))
  return(p)
}
