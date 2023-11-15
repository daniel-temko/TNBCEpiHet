library(ggplot2)
library(ggrepel)
library(cowplot)
library(cluster)
library(foreach)
library(doMC)
registerDoMC(22)
source("../utilities/utils.R")

setwd("../")

# Load data -----------------------------------------------------------------------------------------------------------------------

load(file = "R_Data/Heatmap_Metadata.RData")
load(file = "R_Data/CyTOF_Data.RData")
load(file = "R_Data/CyTOF_Metadata.RData")
#load("R_Data/scRNA_MC_Libraries.RData")
cytof_sample_sheet <- cytof_sample_sheet[which(cytof_sample_sheet$include == "yes"),]

x_shift_res <- read.table("analysis/cytof/x_shift/x_shift_20200609/csv/x_shift_20200609.csv",
                          sep = ",",
                          header = TRUE,
                          stringsAsFactors = FALSE)

# Aggregate clusters -----------------------------------------------------------------------------------------------------------------------

x_shift_res$sample <- factor(x_shift_res$File.Name, 
                             levels = cytof_sample_sheet$gated.prefix,
                             labels = cytof_sample_sheet$Sample)
x_shift_clust <- as.data.frame.matrix(table(x_shift_res$sample, x_shift_res$ClusterID))
x_shift_prop <- t(apply(x_shift_clust, 1, function(x) x/sum(x)))
max_sample1 <- 1000
max_sample2 <- 200

# Calculate diversity -----------------------------------------------------------------------------------------------------------------------

sample_div <- data.frame(sample.name = cytof_sample_sheet$Sample, 
                         cell.line = cytof_sample_sheet$cell.line)

## X-shift shannon diversity
sample_div$xshift.shannon.idx <- apply(x_shift_prop, 1, function(x) {
  -1 * sum(x*log(x), na.rm = TRUE)
})

# Assess Mean Distance (full data)
sample_div$mean.dist <- NA
for(i in 1:nrow(sample_div)){
  sample <- sample_div$sample.name[i]
  print(sample)
  df <- cytof_norm[which(cytof_metadata$sample == sample),]
  av_profile <- apply(df, 2, mean)
  mean_dist <- mean(sqrt(apply(apply(df, 1, function(x) x - av_profile)^2, 2, sum)))
  sample_div$mean.dist[i] <- mean_dist
}

## down-sample
set.seed(20200610)
n_cells <- as.numeric(table(cytof_metadata$sample))
n_ds <- min(min(n_cells), max_sample1)
pw_dist_ids <- lapply(n_cells, function(x) sample(x, n_ds, replace = TRUE))

# Assess Mean Pairwise Distance
#load(file = "R_Data/CyTOF_PW_Dist.RData")
pw_dist_res = foreach(i = 1:nrow(sample_div)) %dopar% {
  sample <- as.character(sample_div$sample.name[i])
  df <- cytof_norm[which(cytof_metadata$sample == sample),]
  df <- df[pw_dist_ids[[i]],]
  pairwise_distances(t(df))
}
names(pw_dist_res) <- as.character(sample_div$sample.name)
save(pw_dist_res, file = "R_Data/CyTOF_PW_Dist.RData")

sample_div$mean.pw.dist.ds <- sapply(pw_dist_res, mean)

# Assess Mean Distance
sample_div$mean.dist.ds <- NA
for(i in 1:nrow(sample_div)){
  sample <- sample_div$sample.name[i]
  print(sample)
  df <- cytof_norm[which(cytof_metadata$sample == sample),]
  df <- df[pw_dist_ids[[i]],]
  
  av_profile <- apply(df, 2, mean)
  mean_dist <- mean(sqrt(apply(apply(df, 1, function(x) x - av_profile)^2, 2, sum)))
  sample_div$mean.dist.ds[i] <- mean_dist
}

## down-sample again
cytof_sil <- cytof_sample_sheet[match(unique(cytof_sample_sheet$cell.line), 
                                             cytof_sample_sheet$cell.line),]

set.seed(20200610)
n_cells <- as.numeric(table(cytof_metadata$sample))
n_ds <- max_sample2
pw_dist_ids <- lapply(n_cells, function(x) sample(x, n_ds, replace = TRUE))
names(pw_dist_ids) <- cytof_sample_sheet$Sample

## Assess silhouette widths
cytof_df <- c()
for(i in 1:nrow(cytof_sil)){
  sample <- cytof_sil$Sample[i]
  print(sample)
  df <- cytof_norm[which(cytof_metadata$sample == sample),]
  df <- df[pw_dist_ids[[sample]],]
  cytof_df <- rbind(cytof_df, df)
}

cytof_d_mat <- dist(cytof_df)
cytof_samp_ids <- rep(1:nrow(cytof_sil), each = max_sample2)
cytof_res_sil <- silhouette(cytof_samp_ids, cytof_d_mat)


# Plot results -----------------------------------------------------------------------------------------------------------------------

pdf("analysis/cytof/diversity/silhouette.pdf")
  plot(cytof_res_sil, 
       col = line_colours$type[cytof_sil$tnbc.subtype],
       main = "Silhouette Widhts")
dev.off()

sil_wid <- data.frame(width = summary(cytof_res_sil)$clus.avg.widths,
                      line = cytof_sil$cell.line)
sil_wid$plot_types <- cytof_sample_sheet$tnbc.subtype[match(sil_wid$line, 
                                                    cytof_sample_sheet$cell.line)]
sil_wid <- sil_wid[order(sil_wid$width),]
pdf("analysis/cytof/diversity/silhouette_barplot.pdf")
  par(mar = c(6, 4, 4, 4))
  barplot(sil_wid$width, names.arg = cytof_sil$cell.line, las = 2, 
          col = line_colours$type[sil_wid$plot_types],
          main = "Average Silhouette Width")
dev.off()

# Save outputs -----------------------------------------------------------------------------------------------------------------------

## Combine replicates
cl_div <- aggregate(cbind(xshift.shannon.idx, mean.dist, mean.dist.ds, mean.pw.dist.ds) ~ cell.line, 
                    data = sample_div, 
                    mean)

cl_div$mean.dist.rank <- get_rank(cl_div$mean.dist)
cl_div$mean.dist.ds.rank <- get_rank(cl_div$mean.dist.ds)
cl_div$mean.pw.dist.ds.rank <- get_rank(cl_div$mean.pw.dist.ds)
cl_div$xshift.shannon.idx.rank <- get_rank(cl_div$xshift.shannon.idx)

cl_div <- cl_div[order(cl_div$mean.pw.dist.ds.rank),]

out <- cl_div[,c("cell.line", 
                 "mean.dist", 
                 "mean.dist.ds", 
                 "mean.pw.dist.ds", 
                 "mean.dist.rank", 
                 "mean.dist.ds.rank", 
                 "mean.pw.dist.ds.rank")]
write.table(out, 
            file = "analysis/cytof/diversity/cytof_diversity.csv",
            sep = ",",
            row.names = FALSE)

out <- cl_div[,c("cell.line",
                 "xshift.shannon.idx",
                 "mean.dist", 
                 "mean.dist.ds", 
                 "mean.pw.dist.ds", 
                 "mean.dist.rank", 
                 "mean.dist.ds.rank", 
                 "mean.pw.dist.ds.rank",
                 "xshift.shannon.idx.rank")]
write.table(out, 
            file = "analysis/cytof/diversity/cytof_diversity_with_xshift.csv",
            sep = ",",
            row.names = FALSE)

save(cl_div, out, file = "R_Data/CyTOF_Diversity.RData")

# Assess correspondence between metrics -----------------------------------------------------------------

measures <- c("mean.dist", "mean.dist.ds", "mean.pw.dist.ds")
plot_list <- list()
plot_num <- 0
for(i in 1:length(measures)){
  for(j in 1:length(measures)){
    plot_num <- plot_num + 1
    print(plot_num)
    Data <- data.frame(cell.line = cl_div$cell.line, 
                       x = cl_div[[measures[j]]], 
                       y = cl_div[[measures[i]]])
    p <- ggplot(Data, aes(x = x, y = y, label = cell.line)) +
        geom_point() + 
        xlab(measures[j]) +
        ylab(measures[i])
    plot_list[[plot_num]] <- p
  }
}
pdf("analysis/cytof/diversity/metric_comparison.pdf")
  plot_grid(plotlist = plot_list)
dev.off()

## Shannon vs. pw dist
pdf("analysis/cytof/diversity/xshift_shannon_pw_dist_ds.pdf", width = 11, height = 8)
  p <- ggplot(cl_div, aes(x = mean.pw.dist.ds, y = xshift.shannon.idx, label = cell.line)) +
    geom_point() + 
    geom_text_repel()
  print(p)
dev.off()
  
## Shannon barplot
plot_df <- cl_div[order(cl_div$xshift.shannon.idx),]
tnbc_types <- as.character(line_annotation[as.character(plot_df$cell.line),c("type")])
pdf("analysis/cytof/diversity/x_shift_shannon_index.pdf", width = 11, height = 8)
  par(mar = c(10,4,4,4))
  barplot(plot_df$xshift.shannon.idx,
          names.arg = plot_df$cell.line,
          col = line_colours$type[tnbc_types], 
          ylab = "xshift.shannon.index",
          las = 2)
dev.off()
