library(ggsci)
library(reshape2)

setwd("../")

out_path <- "analysis/figures/mofa_graphical_abstract/"
cols = pal_d3()(10)[1:6]

# Explicit weights
feature_ids <- list(dataset1 = 1:5, dataset2 = 6:10, dataset3 = 11:15)
n_samples <- 30
n_factors <- 2
n_features <- sum(sapply(feature_ids, length))
n_factor_vals <- n_factors * n_samples
n_feature_vals <- n_features * n_samples

set.seed(1235)
sim_weights <- matrix(rep(0, n_factors * n_features), ncol = n_factors)
f1_nz_ids <- c(1:3, 7, 9, 10)
f2_nz_ids <- c(4, 5, 12, 13, 14)
sim_weights[f1_nz_ids,1] <- rnorm(length(f1_nz_ids))
sim_weights[f2_nz_ids,2] <- rnorm(length(f2_nz_ids))
sim_weights <- pmin(sim_weights, 1.9)
sim_weights <- pmax(sim_weights, -1.9)

sim_noise <- matrix(rnorm(n_feature_vals, sd=0.8), nrow = n_samples)
sim_factor_vals <- matrix(rnorm(n_factor_vals), nrow = n_samples)
sim_exp1 <- sim_factor_vals[,1,drop=FALSE] %*% t(sim_weights[,1]) 
sim_exp2 <- sim_factor_vals[,2,drop=FALSE] %*% t(sim_weights[,2]) 
sim_samples <- sim_factor_vals %*% t(sim_weights) + sim_noise
sim_cov <- pmax((sim_factor_vals[,1] + rnorm(n_samples) + 5) * 5, 1)

r2_1_1 <- 1 - sum((sim_samples - sim_exp1)[,feature_ids[[1]]]^2) / sum(sim_samples[,feature_ids[[1]]]^2)
r2_1_2 <- 1 - sum((sim_samples - sim_exp1)[,feature_ids[[2]]]^2) / sum(sim_samples[,feature_ids[[2]]]^2)
r2_1_3 <- 1 - sum((sim_samples - sim_exp1)[,feature_ids[[3]]]^2) / sum(sim_samples[,feature_ids[[3]]]^2)

r2_2_1 <- 1 - sum((sim_samples - sim_exp2)[,feature_ids[[1]]]^2) / sum(sim_samples[,feature_ids[[1]]]^2)
r2_2_2 <- 1 - sum((sim_samples - sim_exp2)[,feature_ids[[2]]]^2) / sum(sim_samples[,feature_ids[[2]]]^2)
r2_2_3 <- 1 - sum((sim_samples - sim_exp2)[,feature_ids[[3]]]^2) / sum(sim_samples[,feature_ids[[3]]]^2)

r2_sim <- data.frame(factor = rep(c(1, 2), each = 3),
                     dataset = factor(rep(1:3, 2)),
                     r2 = c(r2_1_1, r2_1_2, r2_1_3, r2_2_1, r2_2_2, r2_2_3))

plot_ids <- c(1, 2, n_samples)
s1_vals <- sim_samples[plot_ids[1],]
s2_vals <- sim_samples[plot_ids[2],]
s3_vals <- sim_samples[plot_ids[3],]

fn <- paste0("f1.pdf")
pdf(paste0(out_path, fn), width = 8, height = 3)
par(mfrow = c(1, 3), mar = c(8, 6, 4, 4))
for(i in 1:length(feature_ids)){
  x <- barplot(sim_weights[feature_ids[[i]],1], 
               col = cols[i], 
               ylim = c(-2, 2),
               ylab = "Factor 1 Loading", 
               main = paste0("Dataset ", i),
               cex.lab = 1.5,
               cex.main = 1.5)
  axis(side = 1, 
       labels = paste0("gene", 1:length(feature_ids[[i]])), 
       at = x[,1], 
       las = 2, 
       cex.axis = 1.5)
}
dev.off()

fn <- paste0("f2.pdf")
pdf(paste0(out_path, fn), width = 8, height = 3)
par(mfrow = c(1, 3), mar = c(8, 6, 4, 4))
for(i in 1:length(feature_ids)){
  x <- barplot(sim_weights[feature_ids[[i]],2], 
               col = cols[i], 
               ylim = c(-2, 2), 
               ylab = "Factor 2 Loading", 
               main = paste0("Dataset ", i),
               cex.lab = 1.5,
               cex.main = 1.5)
  axis(side = 1, 
       labels = paste0("gene", 1:length(feature_ids[[i]])), 
       at = x[,1], 
       las = 2, 
       cex.axis = 1.5)
}
dev.off()

fn <- paste0("s1_vals.pdf")
pdf(paste0(out_path, fn), width = 8, height = 3)
par(mfrow = c(1, 3), mar = c(8, 6, 4, 4))
for(i in 1:length(feature_ids)){
  barplot(t(s1_vals)[, feature_ids[[i]]], 
          col = cols[i],
          ylim = c(-5, 5),
          ylab = paste0("Sample ", plot_ids[1], " Change vs. Average"),
          main = paste0("Dataset ", i),
          cex.lab = 1.5,
          cex.main = 1.5)
  axis(side = 1, 
       labels = paste0("gene", 1:length(feature_ids[[i]])), 
       at = x[,1], 
       las = 2, 
       cex.axis = 1.5)
}
dev.off()

fn <- paste0("s2_vals.pdf")
pdf(paste0(out_path, fn), width = 8, height = 3)
par(mfrow = c(1, 3), mar = c(8, 6, 4, 4))
for(i in 1:length(feature_ids)){
  barplot(t(s2_vals)[, feature_ids[[i]]], 
          col = cols[i],
          ylim = c(-5, 5),
          ylab = paste0("Sample ", plot_ids[2], " Change vs. Average"),
          main = paste0("Dataset ", i),
          cex.lab = 1.5,
          cex.main = 1.5)
  axis(side = 1, 
       labels = paste0("gene", 1:length(feature_ids[[i]])), 
       at = x[,1], 
       las = 2, 
       cex.axis = 1.5)
}
dev.off()

fn <- paste0("s3_vals.pdf")
pdf(paste0(out_path, fn), width = , height = 3)
par(mfrow = c(1, 3), mar = c(8, 6, 4, 4))
for(i in 1:length(feature_ids)){
  barplot(t(s3_vals)[, feature_ids[[i]]], 
          col = cols[i],
          ylim = c(-5, 5),
          ylab = paste0("Sample ", plot_ids[3], " Change vs. Average"),
          main = paste0("Dataset ", i),
          cex.lab = 1.5,
          cex.main = 1.5)
  axis(side = 1, 
       labels = paste0("gene", 1:length(feature_ids[[i]])), 
       at = x[,1], 
       las = 2, 
       cex.axis = 1.5)
}
dev.off()

sim_scores <- as.data.frame(sim_factor_vals)
colnames(sim_scores) <- c("factor.1", "factor.2")
sim_scores$sample <- c("sample.1", rep("", n_samples - 2), paste0("sample.", n_samples))

fn <- paste0("sim_scores.pdf")
pdf(paste0(out_path, fn))
ggplot(sim_scores, aes(x = factor.1, y = factor.2, color = sample)) +
  geom_point() + 
  scale_color_manual(values = c("#000000", cols[4:5])) +
  theme(text = element_text(size = 30), legend.position = "none") +
  geom_text_repel(aes(label = sample), size = 8) +
  xlab("Factor 1") +
  ylab("Factor 2")
dev.off()

sub <- subset(r2_sim, factor == 1)
fn <- paste0("f1_var_exp.pdf")
pdf(paste0(out_path, fn), width = 4, height = 2.5)
ggplot(sub, aes(x = dataset, y = r2, fill = dataset)) +
  geom_bar(stat = "identity") + 
  ylim(c(0,1)) + 
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values = cols[1:3]) +
  ggtitle("Var. Exp. Factor 1 ") + 
  ylab("Var. Explained")
dev.off()

sub <- subset(r2_sim, factor == 2)
fn <- paste0("f2_var_exp.pdf")
pdf(paste0(out_path, fn), width = 4, height = 2.5)
ggplot(sub, aes(x = dataset, y = r2, fill = dataset)) +
  geom_bar(stat = "identity") + 
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values = cols[1:3]) + 
  ylim(c(0,1)) +
  ggtitle("Var. Exp. Factor 2") + 
  ylab("Var. Explained")
dev.off()

df <- data.frame(factor.1 = sim_factor_vals[,1], cov = sim_cov) 
fn <- paste0("f1_covariate.pdf")
pdf(paste0(out_path, fn), width = 4, height = 2.5)
ggplot(df, aes(x = factor.1, y=cov)) +
  geom_point() +
  xlab("Factor 1") +
  ylab("Covariate")
dev.off()