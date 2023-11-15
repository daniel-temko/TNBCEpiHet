library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(cowplot)

#' Generate colour gradient for plotting based on a give base color
GenerateColourRamp <- function(base_col, bias = 1){
  mid_col <- colorRampPalette(colors = c("#FFFFFF", base_col), space = "Lab")(100)[50]
  ColourRamp <- colorRampPalette(colors = c("#FFFFFF", mid_col, base_col), 
                                 space = "Lab", 
                                 bias = bias)
  return(ColourRamp(100))
}

#' Plot a data heatmap. Mean-centers the data as a pre-processing step. Clusters using Euclidean 
#' distance and Ward.D2 for columns and Euclidean distance or correlation distance -- 
#' Euclidean distance after scaling rows -- with Ward.D2 for rows. If cap is TRUE, then
#' values are capped at +/-3 for visualization
PlotHeatmap <- function(data, annotation_col, annotation_colours,
                        method = "euclidean",
                        colour = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                        cluster_rows = TRUE,
                        cap = TRUE,
                        label_samples = TRUE,
                        cellwidth = NA,
                        cellheight = NA){
  
  # center (for visualisation, and for correlation distance, if applicable)
  data <- t(apply(data, 1, function(x) x - mean(x)))
  
  if (method == "euclidean") {
    col_dist <- dist(t(data), method = "euclidean")  
  } else if (method == "correlation") {
    data_norm <- scale(data)
    col_dist <- dist(t(data_norm), method = "euclidean")
  } else {
    stop("incorrect method specification")
  }
  col_hclust <- hclust(col_dist, method = "ward.D2")
  
  if(cluster_rows){
    row_dist <- dist(data, method = "euclidean")
    row_hclust <- hclust(row_dist, method = "ward.D2")
  } else {
    row_hclust <- FALSE
  }
  
  if (cap) {
    cap_val <- 3
    data <- t(apply(data, 1, function(x) sapply(x , function(y) min(max(-cap_val,y),cap_val))))
  }
  
  if(nrow(data) <= 100) {
    show_rownames <- TRUE
  } else {
    show_rownames <- FALSE
  }
  
  if(ncol(data) >= 40) {
    fontsize_col <- 400 / ncol(data)
  } else {
    fontsize_col <- 10
  }
  
  pheatmap(data, 
           annotation_col = annotation_col, 
           annotation_colors = annotation_colours,
           color = colour,
           cluster_cols = col_hclust,
           cluster_rows = row_hclust,
           show_rownames = show_rownames,
           fontsize_col = fontsize_col,
           treeheight_row = 0,
           treeheight_col = 100,
           show_colnames = label_samples,
           cellwidth = cellwidth,
           cellheight = cellheight)   
  
  return(col_hclust)
}

#' Plots raster of scatter plots of samples ordinated by first three
#' principal components. Samples are coloured using colours argument
PlotPCARaster <- function(data, colours, file_name){
  raster_plots <- list()
  counter <- 1
  for(j in 1:3){
    for(i in 1:3) {
      if(i != j) {
        plot_data <- data.frame(x = data[[paste0("PC", i)]], y = data[[paste0("PC", j)]])
        p <- ggplot(data = plot_data, aes(x = x, y = y)) + 
          geom_point(color = colours) +
          xlab(paste0("PC", i)) +
          ylab(paste0("PC", j))
        raster_plots[[counter]] <- p  
      }
      counter <- counter + 1
    }
  }
  pdf(file = file_name)
  p <- do.call(plot_grid, c(raster_plots, ncol = 3))
  print(p)
  dev.off()
}

#' Wrapper function that analyses a range of data cuts with different 
#' numbers of HVG's, plots heatmaps, runs PCA, and plots the results
GenerateExploratoryPlots <- function(data, 
                                     top_features, 
                                     n_pcs, 
                                     base_col,
                                     base_dir,
                                     perc_features = c(seq(10, 50, 10), 100),
                                     annotation_col = line_annotation,
                                     annotation_colours = line_colours,
                                     heatmap_col = NULL,
                                     heatmap_colours = NULL,
                                     umap = FALSE,
                                     pca = TRUE,
                                     max_features = 25000,
                                     label_samples = TRUE,
                                     plot_pdf = TRUE,
                                     cap = 3){
  if(is.null(heatmap_col)) heatmap_col <- annotation_col
  if(is.null(heatmap_colours)) heatmap_colours <- annotation_colours
  
  if(length(perc_features) != length(n_pcs)) {
    warning("wrong size n_pcs")
    return(1)
  }
  
  orig_dir <- getwd()
  if(!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)
  setwd(base_dir)
  
  if(!dir.exists(paste0("pca_data"))) dir.create(paste0("pca_data"), recursive = TRUE)
  if(!dir.exists(paste0("pca_summary"))) dir.create(paste0("pca_summary"), recursive = TRUE)
  if(!dir.exists(paste0("raw_data"))) dir.create(paste0("raw_data"), recursive = TRUE)
  if(!dir.exists(paste0("r_data"))) dir.create(paste0("r_data"), recursive = TRUE)
  if(plot_pdf){
    if(!dir.exists(paste0("raw_data/pdf"))) dir.create(paste0("raw_data/pdf"), recursive = TRUE)
  }
  
  colour_ramp <- GenerateColourRamp(base_col)
  
  data <- na.omit(data)
  
  types <- as.character(annotation_col$type[match(colnames(data), rownames(annotation_col))])
  
  # calculate feature numbers for each pass, floored at two
  num_features <- round(perc_features * (1/100) * length(top_features))
  
  pretty_perc <- sapply(perc_features, function(x) formatC(x, width = 2, flag = 0))
  
  # update top feature list to conform with row names
  top_features <- c(top_features[which(top_features %in% rownames(data))],
                    rownames(data)[which(!(rownames(data) %in% top_features))])
  
  for(i in seq_along(perc_features)){
    print(i)
    
    if(num_features[i] < 3) {next}
    
    data_top_features <- data[top_features[1:num_features[i]],] 
    
    # Raw heatmaps (skip if num_features > 15000)
    if(num_features[i] <= max_features){
      heatmap_file_name <- paste0("raw_data/heatmap_euc_top", pretty_perc[i], "pc.png")
      clustering_file_name <- paste0("r_data/raw_data_euc_top", pretty_perc[i], "pc_clustering.RData")
      png(file = heatmap_file_name, res = 300, width = 9, height = 9, units = "in", type = "cairo")
      par(mar = c(4,4,4,4))
      col_hclust <- PlotHeatmap(data_top_features, 
                                heatmap_col, 
                                heatmap_colours, 
                                colour = colour_ramp,
                                label_samples = label_samples,
                                cap = cap)
      dev.off()
      save(col_hclust, file = clustering_file_name)
      
      if(plot_pdf){
        heatmap_file_name <- paste0("raw_data/pdf/heatmap_euc_top", pretty_perc[i], "pc.pdf")
        pdf(heatmap_file_name)
        PlotHeatmap(data_top_features, 
                    heatmap_col, 
                    heatmap_colours, 
                    colour = colour_ramp,
                    label_samples = label_samples,
                    cap = cap)
        dev.off()
      }
      
      heatmap_file_name <- paste0("raw_data/heatmap_cor_top", pretty_perc[i], "pc.png")
      clustering_file_name <- paste0("r_data/raw_data_cor_top", pretty_perc[i], "pc_clustering.RData")
      png(file = heatmap_file_name, res = 300, width = 9, height = 9, units = "in", type = "cairo")
      par(mar = c(4,4,4,4))
      col_hclust <- PlotHeatmap(data_top_features, 
                                heatmap_col, 
                                heatmap_colours, 
                                colour = colour_ramp,
                                method = "correlation",
                                label_samples = label_samples,
                                cap = cap)
      dev.off()
      save(col_hclust, file = clustering_file_name)
      
      if(plot_pdf){
        heatmap_file_name <- paste0("raw_data/pdf/heatmap_cor_top", pretty_perc[i], "pc.pdf")
        pdf(heatmap_file_name)
        PlotHeatmap(data_top_features, 
                    heatmap_col, 
                    heatmap_colours, 
                    colour = colour_ramp,
                    method = "correlation",
                    label_samples = label_samples,
                    cap = cap)
        dev.off()
      }
    }
    
    #if(num_features[i] > 100000) {next}
    
    if(umap){
      # Raw umap
      print('UMAP...')
      set.seed(1234)
      data_umap <- umap(t(data_top_features), n_neighbors = 5)$layout
      colnames(data_umap) <- c("x", "y")
      
      pdf(file = paste0("raw_data/umap_top", pretty_perc[i], "pc.pdf"))
      p <- ggplot(data = as.data.frame(data_umap), aes(x = x, y = y)) + 
        geom_point(color = annotation_colours$type[types])
      if(label_samples){
        p <- p + geom_text_repel(label = rownames(data_umap), max.overlaps = Inf)
      }
      print(p)
      dev.off()
    }
    #if(num_features[i] > 50000) {next}
    
    if(pca){
      print("PCA...")
      data_pca <- prcomp(t(as.matrix(data_top_features)))
      
      # PC1 vs PC2 scatter
      pdf(file = paste0("pca_data/pca1-pca2-top", pretty_perc[i], "pc.pdf"))
      p <- ggplot(data = as.data.frame(data_pca$x), aes(x = PC1, y = PC2)) + 
        #geom_point(color = annotation_colours$type[types])
        geom_point(aes(color = types)) +
        scale_color_manual(values = annotation_colours$type)
      if(label_samples){
        p <- p + geom_text_repel(label = rownames(data_pca$x), max.overlaps = Inf)
      }
      print(p)
      dev.off()
      
      # PC1-3 raster
      PlotPCARaster(as.data.frame(data_pca$x), 
                    annotation_colours$type[types], 
                    paste0("pca_data/pca1to3_top", pretty_perc[i], ".pdf"))
      
      # PCA variance explained
      pdf(file = paste0("pca_summary/pca_variance_explained_top", pretty_perc[i], "pc.pdf"))
      plot(summary(data_pca)$importance[2,], 
           xlab = "Component", 
           ylab = "Variance Explained")
      dev.off()
      
      data_pca_top_dims <- t(data_pca$x[,1:min(n_pcs[i], ncol(data_pca$x))])
      
      # PCA heatmaps
      heatmap_file_name <- paste0("pca_data/heatmap_top", pretty_perc[i], "pc_ndim", n_pcs[i], ".pdf")
      clustering_file_name <- paste0("r_data/pca_data_top", pretty_perc[i], "pc_clustering.RData")
      pdf(heatmap_file_name)
      col_hclust <- PlotHeatmap(data_pca_top_dims, 
                                heatmap_col, 
                                heatmap_colours, 
                                cluster_rows = FALSE, 
                                cap = FALSE,
                                label_samples = label_samples)
      dev.off()
      save(col_hclust, file = clustering_file_name)
    }
    
    # PCA umap
    if(umap){
      set.seed(1234)
      data_pca_umap <- umap(t(data_pca_top_dims), n_neighbors = 5)$layout
      colnames(data_pca_umap) <- c("x", "y")
      
      pdf(file = paste0("pca_data/umap_top_", pretty_perc[i], "pc_ndim", n_pcs[i], ".pdf"))
      p <- ggplot(data = as.data.frame(data_pca_umap), aes(x = x, y = y)) + 
        geom_point(color = annotation_colours$type[types]) 
      if(label_samples){
        p <- p + geom_text_repel(label = rownames(data_pca_umap), max.overlaps = Inf)  
      }
      print(p)
      dev.off()
    }
  }
  setwd(orig_dir)
}