source("../utilities/utils.R")
source("../utilities/sampling.R")
source("../utilities/signatures.R")
source("../utilities/mofa2_factor_calls.R")
source("../src/hexplot/hexplot.R")
library(Seurat)
library(reshape2)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(SingleCellExperiment)
library(foreach)
library(doMC)
registerDoMC(22)

setwd("../")

# Load data ---------------------------------------------------------------------------------------------------------------------------------------

load(file = "R_Data/SC_BC_DR.RData")
load(file = "R_Data/SC_Metadata.RData")
load(file = "R_Data/Heatmap_Metadata.RData")
load(file = "R_Data/Gene_Sigs.RData")
load(file = "R_Data/RNA_Data.RData")
load(file = "R_Data/RNA_Metadata.RData")
load(file = "R_Data/MOFA2_TopN_Heatmap_Feature_Filtering_Seed20200225_RNA_Weights.RData")

# MOFA diversity ---------------------------------------------------------------------------------------------------------------------------------------

# Project onto MOFA tnbc type-linked factors
if(file.exists("R_Data/SC_MOFA_Scores.RData")){
  load("R_Data/SC_MOFA_Scores.RData")
} else{
  keep <- rownames(sc) %in% rownames(rna_weights)
  data_list <- list(rna = sc[["lnbc"]]@data[keep, ])
  weight_list <- list(rna = rna_weights)
  call_data <- prepare_mofa_calls(data_list, weight_list)
  mofa_scores <- call_mofa_factor_scores(call_data)
  orig_idents <- sc$orig.ident
  save(orig_idents, mofa_scores, file = "R_Data/SC_MOFA_Scores.RData")
}

# MOFA ITH analysis
cl <- factor(sc@meta.data$orig.ident)
wss <- apply(mofa_scores, 2, function(x){
  sum(tapply(x, cl, var) * (table(cl) - 1))
})
tss <- apply(mofa_scores, 2, var) * (nrow(mofa_scores) - 1)
mofa_het <- data.frame(wss = wss, tss = tss)
mofa_het$score <- mofa_het$wss / mofa_het$tss
mofa_het$factor <- colnames(mofa_scores)

pp <- mofa_het
pp$factor <- factor(pp$factor, levels = pp$factor[order(pp$score)])
pdf("analysis/scrna_seq/mofa_projection/within_ss.pdf")
ggplot(pp, aes(x = factor, y = score, fill = factor)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  xlab("Factor") +
  ylab("Intra-cell-line variance proportion") + 
  my_theme()
dev.off()
