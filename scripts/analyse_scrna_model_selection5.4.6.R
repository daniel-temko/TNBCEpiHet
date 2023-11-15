source("../utilities/utils.R")
library(mclust)
library(ggplot2)
library(Matrix)
library(reshape2)
library(ggrepel)

setwd("../")

# Load data ---------------------------------------------------------------------------------------------------------------------------------------

obj.sc <- load(file = "R_Data/SC_Metadata.RData")
obj.mo <- load(file = "R_Data/SC_MOFA_Scores.RData")
obj.s3 <- load(file = "R_Data/SC_Sample3.RData")

# Clustering ---------------------------------------------------------------------------------------------------------------------------------------

# Assess evidence for MOFA score clusters
if(file.exists("R_Data/SC_MOFA_LL_5.4.6.RData")){
  load("R_Data/SC_MOFA_LL_5.4.6.RData")
} else {
  ll <- list()
  for(i in 1:nrow(scrna_sample_table)){
    sample <- as.character(scrna_sample_table$sample3[i])
    print(sample)
    mm <- mofa_scores[sample3 == sample, c(2, 3, 6)]
    ll[[sample]] <- sapply(1:5, function(x) Mclust(mm, x, modelNames = "VVV")$loglik)
  }
  save(ll, file = "R_Data/SC_MOFA_LL_5.4.6.RData")  
}

ll2 <- sapply(ll, function(x) x[2:5]-x[1:4])
p <- melt(cbind(as.data.frame(ll2), k = 2:5), id.var = "k")
colnames(p) <- c("Clusters", "Sample", "DeltaLogLik") 
p$Sample <- factor(p$Sample, levels = sort(scrna_sample_table$sample3))
p$Label <- as.character(p$Sample)
p$Label[which(p$Clusters != 2 | !p$Sample %in% c("HDQP1", "SUM185", "EMG3"))] = ""
pdf("analysis/scrna_seq/mofa_projection/clustering/delta_loglik_vs_last_5.4.6.pdf")
ggplot(p, aes(x = Clusters, y = DeltaLogLik, color = Sample)) + 
  geom_point() +
  geom_text_repel(aes(label = Label), max.overlaps = Inf, show.legend = FALSE) +
  ylab("Change in log likelihood") +
  my_theme() +
  theme(legend.position="top", legend.title = element_blank(), legend.text = element_text(size = 10))
dev.off()