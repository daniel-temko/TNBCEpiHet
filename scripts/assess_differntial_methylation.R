library(limma)
library(dmrc)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

setwd("../")

load("R_Data/Methyl_Regions.RData")
load("R_Data/Methyl_Data.RData")

ann_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

type <- "mesenchymal"

tnbc_type <- factor(sapply(methyl_metadata$tnbc.subtype, function(x){
  if (x == type){
    type
  } else {
    "other"
  }
}), levels = c(type, "other"))

#tnbc_type <- factor(methyl_metadata$tnbc.subtype, 
#                    levels = sort(unique(methyl_metadata$tnbc.subtype)))

design <- model.matrix(~0 + tnbc_type, data = methyl_metadata)
colnames(design) <- levels(tnbc_type)

fit <- lmFit(m_vals, design)

cont_matrix <- makeContrasts(paste0(type, " - other"), levels = design)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))

all(rownames(m_vals) %in% ann_450k$Name)
ann_450k_sub <- ann_450k[match(rownames(m_vals), ann_450k$Name), c(1:4,12:19,24:ncol(ann_450k))]

dmps <- topTable(fit2, num=Inf, coef=1, genelist = ann_450k_sub)

#' assumes region.id is a factor
aggregate_dmp_df <- function(dmp_df, map, id){
  map <- map[which(map$locus %in% rownames(dmp_df)),]
  locus_vals <- dmp_df[as.character(map$locus), ]
  
  !("region.id" %in% colnames(locus_vals))
  locus_vals$region.id <- map$region.id
  cts <- aggregate(adj.P.Val ~ region.id, locus_vals, length, drop = FALSE)[,2]
  cts[which(is.na(cts))] <- 0
  min_padj <- aggregate(adj.P.Val ~ region.id, locus_vals, min, drop = FALSE)[,2]
  res <- data.frame(row.names = levels(map$region.id), n = cts, min_padj = min_padj)
  colnames(res) <- paste0(colnames(res), "_", id)
  return(res)
}


# Summarize to TSS, GB, and SE levels
summarize_dmps_to_regions <- function(dmps, map, alpha = 0.05){
  map <- map[which(map$locus %in% rownames(dmps)),]
  map$region.id <- factor(map$region.id)
  
  sig_up <- subset(dmps, (adj.P.Val < alpha) & (logFC > 0))
  sig_dn <- subset(dmps, (adj.P.Val < alpha) & (logFC < 0))
  no_sig <- subset(dmps, (adj.P.Val) >= alpha)
  
  sig_up_agg <- aggregate_dmp_df(sig_up, map, id = "sig_up")
  sig_dn_agg <- aggregate_dmp_df(sig_dn, map, id = "sig_dn")
  no_sig_agg <- aggregate_dmp_df(no_sig, map, id = "non_sig")
  res <- cbind(sig_up_agg, sig_dn_agg, no_sig_agg)
  
  res$n_tot <- as.numeric(table(map$region.id))
  stopifnot(all(res$n_sig_up + res$n_sig_dn + res$n_non_sig == res$n_tot))
  
  res$prop_sig_up <- res$n_sig_up / res$n_tot
  res$prop_sig_dn <- res$n_sig_dn / res$n_tot
  
  locus_vals <- dmps[as.character(map$locus), ]
  !("region.id" %in% colnames(locus_vals))
  locus_vals$region.id <- map$region.id
  res$av_lfc <- aggregate(logFC ~ region.id, locus_vals, mean, drop = FALSE)[,2]
  res <- res[,c("av_lfc", "n_sig_up", "n_sig_dn", "n_non_sig", "n_tot", "prop_sig_up", 
                "prop_sig_dn", "min_padj_sig_up", "min_padj_sig_dn")]
  res <- res[order(pmin(res$min_padj_sig_dn, res$min_padj_sig_up)), ]
  return(res)
}


