#####
# analyze each sample using Negative Binomial distribution
#####



################################################################################
# load library
################################################################################
library(Seurat)
library(ggplot2)

################################################################################
# load data
################################################################################
obj_phe <- load("R_Data/SC_Metadata.RData")
obj_met <- load("R_Data/Heatmap_Metadata.RData")
obj_dat <- load("R_Data/SC_Filt.RData")

hkg <- c("ACTB", "GAPDH", "C1orf43", "CHMP2A", "GPI", "PSMB2", "PSMB4", "RAB7A",
         "REEP5", "SNRPD3", "VCP", "VPS29")

scale_fac <- 10000

################################################################################
# transform data
################################################################################

exp_num <- rowSums(sc[["RNA"]]@counts > 0)
keep <- exp_num >= 10
sc <- sc[keep, ]

################################################################################
# local functions
################################################################################

est_var <- function(exp_count_mat, scale_factor = scale_fac){
  size_factors <- colSums(exp_count_mat) / scale_factor
  exp_norm_mat <- sweep(exp_count_mat, 2, size_factors, "/")
  gene_cs_vars <- apply(exp_norm_mat, 1, var)
  gene_cs_means <- rowMeans(exp_norm_mat)
  bias_factors <- (gene_cs_means / ncol(exp_count_mat)) * sum(1 / size_factors)
  est_raw_var <- gene_cs_vars - bias_factors
  list(cells = data.frame(size_factors = size_factors),
       genes = data.frame(cs_mean = gene_cs_means, 
                          cs_var = gene_cs_vars, 
                          est_raw_var = est_raw_var))
}

################################################################################
# exploratory
################################################################################

pdf("analysis/scrna_seq/diversity/common_scale_variance.pdf")
par(mfrow = c(2, 2))
lapply(levels(scrna_sample_table$sample3), function(x){
  message(x)
  ids <- sc$sample3 == x
  x_type <- scrna_sample_table$tnbc.subtype[match(x, scrna_sample_table$sample3)]
  sub <- sc[, ids]
  exp_count_mat <- sub[["RNA"]]@counts
  size_factors <- colSums(exp_count_mat) / scale_fac
  exp_norm_mat <- sweep(exp_count_mat, 2, size_factors, "/")
  gmeans <- rowMeans(exp_norm_mat)
  gvars <- apply(exp_norm_mat, 1, var)
  pois_vars <- gmeans / mean(size_factors)
  plot(log(gmeans), log(gvars), cex = 0.1, xlim = c(-10, 5), ylim = c(-10, 5),
       main = x)
  lines(log(gmeans)[order(log(gmeans))], 
        log(pois_vars)[order(log(gmeans))], 
        col = line_colours$type[x_type])
})
dev.off()

################################################################################
# analysis of variance
################################################################################

## raw variance estimates
if(file.exists("R_Data/SC_Var_Est.RData")){
  load("R_Data/SC_Var_Est.RData")
} else {
  var_est <- lapply(levels(scrna_sample_table$sample3), function(x){
    message(x)
    ids <- sc$sample3 == x
    sub <- sc[, ids]
    est_var(sub[["RNA"]]@counts)
  })
  names(var_est) <- levels(scrna_sample_table$sample3)
  save(var_est, file = "R_Data/SC_Var_Est.RData")
}

## infer theta vals using glm approach
pdf("analysis/scrna_seq/diversity/common_scale_variance_fit.pdf")
par(mfrow = c(2, 2))
theta_est <- lapply(names(var_est), function(x) {
  message(x)
  x_type <- scrna_sample_table$tnbc.subtype[match(x, scrna_sample_table$sample3)]
  vv <- var_est[[x]]
  mean_s <- mean(vv$cells$size_factors)
  cs_means <- vv$genes$cs_mean + eps
  cs_vars <- vv$genes$cs_var + eps
  offs <- cs_means / mean_s
  x_vals <- (cs_means) ^ 2
  #res <- glm(cs_vars ~ x_vals + 0, offset = offs, family = Gamma(link = "identity"))
  res <- lm(cs_vars ~ x_vals + 0, offset = offs)
  theta_hat <- res$coefficients[1]
  #plot(log(cs_means), log(cs_vars), cex = 0.1, pch = 16, xlim = c(-10, 5), ylim = c(-10, 5),
  #     main = x)
  plot(log(cs_means), log(cs_vars), cex = 0.1, pch = 16,
       main = x)
  abline(a = -log(mean_s), b = 1, col = line_colours$type[x_type])
  abline(a = log(theta_hat), b = 2, col = line_colours$type[x_type])
  theta_hat
})
dev.off()
names(theta_est) <- names(var_est)

pdf("analysis/scrna_seq/diversity/theta_comp.pdf")
gg <- data.frame(sample = names(theta_est),
                 theta_est = sapply(theta_est, function(x) x))
gg$tnbc.type <- scrna_sample_table$tnbc.subtype[match(gg$sample, scrna_sample_table$sample3)]
gg$sample <- factor(gg$sample, levels = gg$sample[order(gg$theta_est)])
ggplot(gg, aes(x = sample, y = theta_est, fill = tnbc.type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = line_colours$type) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
kruskal.test(theta_est ~ tnbc.type, gg)
dev.off()
## transform variance estimates to dispersion estimates (alpha) and plot
pdf("analysis/scrna_seq/diversity/dispersion_estimates_cpht1.pdf")
par(mfrow = c(2, 2))
disp_est <- lapply(names(var_est), function(x){
  message(x)
  x_type <- scrna_sample_table$tnbc.subtype[match(x, scrna_sample_table$sample3)]
  vv <- var_est[[x]]
  eps <- 1e-9
  dd <- pmax(vv$est_raw_var / vv$cs_mean^2, eps)
  keep <- vv$cs_mean > 1 & dd < 1
  plot(log(vv$cs_mean)[keep], dd[keep], ylim = c(0, 1))
  #comm_disp <- mean(disp_est[keep])
  res <- glm(dd[keep] ~ 1, family = Gamma(link = "identity"))
  abline(h = res$coefficients[1], col = line_colours$type[x_type])
  res$coefficients[1]
})
dev.off()
names(disp_est) <- names(var_est)

gg <- data.frame(row.names = names(disp_est),
                 sample = names(disp_est),
                 disp = sapply(disp_est, function(x) x),
                 type = scrna_sample_table$tnbc.subtype[match(names(disp_est), scrna_sample_table$sample3)])
gg$sample <- factor(gg$sample, levels = gg$sample[order(gg$disp)])
ggplot(gg, aes(x = sample, y = disp, fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = line_colours$type) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
kruskal.test(disp ~ type, gg)

## inference of common theta
a <- var_est$SUM1315
mean_s <- mean(a$cells$size_factors)
cs_means <- a$genes$cs_mean + eps
cs_vars <- a$genes$cs_var + eps
offs <- cs_means / mean_s + eps
x_vals <- (cs_means)^2
y_vals <- cs_vars + eps
res <- glm(y_vals ~ x_vals + 0,
           offset = offs, 
           family = Gamma(link = "identity"))
plot(log(a$genes$cs_mean), log(a$genes$cs_var), cex = 0.1, pch = 16)
abline(a = -log(mean_s), b = 1)
abline(a = log(res$coefficients[1]), b = 2)

## transform variance estimates to dispersion estimates (alpha) and plot
eps <- 1e-5
a <- var_est$SUM1315
x_vals <- 1 / (a$cs_mean)
y_vals <- pmax(a$est_raw_var / a$cs_mean^2, eps)
fit <- glm(y_vals ~ x_vals, family = Gamma(link = "identity"))
pred_vals <- predict(fit, data.frame(x_vals = x_vals))
plot(log(a$cs_mean), log(a$est_raw_var / a$cs_mean^2))
lines(log(a$cs_mean)[order(log(a$cs_mean))], log(pred_vals)[order(log(a$cs_mean))], 
      col = "blue")
#lines(log(log(a$cs_mean)[order(log(a$cs_mean))]), log(pred_vals)[order(log(a$cs_mean))], col = "blue")
plot(log(a$cs_mean), a$est_raw_var / a$cs_mean^2, ylim = c(0, 1))

################################################################################
# house-keeper gene analysis
################################################################################

## combine raw variance for housekeeping genes into matrix and analyze
hkg_pres <- hkg[which(hkg %in% rownames(sc))]
hkg_var <- sapply(var_est, function(x) x$gene[hkg_pres, "est_raw_var"])
rownames(hkg_var) <- hkg_pres
library(reshape)
hkg_df <- melt(hkg_var)
hkg_df$type <- scrna_sample_table$tnbc.subtype[match(hkg_df[,2], scrna_sample_table$sample3)]

ggplot(hkg_df, aes(x = type, y = log(value), fill = type)) +
  geom_boxplot() +
  scale_fill_manual(values = line_colours$type) +
  facet_grid(. ~ X1)

hkg_mean <- sapply(var_est, function(x) x$gene[hkg_pres, "cs_mean"])
rownames(hkg_mean) <- hkg_pres
library(reshape)
hkg_df <- melt(hkg_mean)
hkg_df$type <- scrna_sample_table$tnbc.subtype[match(hkg_df[,2], scrna_sample_table$sample3)]

ggplot(hkg_df, aes(x = type, y = log(value), fill = type)) +
  geom_boxplot() +
  scale_fill_manual(values = line_colours$type) +
  facet_grid(. ~ X1)
