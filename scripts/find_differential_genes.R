library(DESeq2)

setwd("../")

load("R_Data/RNA_Data.RData")
load("R_Data/TF_Data.RData")
load("R_Data/TF_UT_Data.RData")

design <- data.frame(group = c("mesenchymal", "basal", "luminal"))

#########################################################################################################
# Differentially expressed genes
#########################################################################################################

if(!dir.exists("analysis/rna_seq/differential_genes")) {dir.create("analysis/rna_seq/differential_genes")}
for(i in 1:nrow(design)){
  print(i)
  
  condition <- factor(sapply(rna_metadata$tnbc.subtype, function(x) if(x == design$group[i]) { x } else { "other" }))
  sample.sheet <- data.frame(sample = colnames(rna_counts), condition = condition)
  
  dds <- DESeqDataSetFromMatrix(countData =  rna_counts, 
                                colData = sample.sheet, 
                                design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > 0, ]

  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", as.character(design$group[i]), "other"))
  
  res$gene.symbol <- rownames(res)
  res$dbd.tf <- res$gene.symbol %in% tf_df$Gene.name
  res$ut.tf <- res$gene.symbol %in% tf_ut$HGNC.symbol
  res <- res[,c((ncol(res) - 2):ncol(res), 1:(ncol(res) - 3))]
  res <- res[order(res$padj),]
  
  sub <- subset(res, (padj < 0.05) & (abs(log2FoldChange) > 1))
  
  assign(paste0(as.character(design$group[i]), "_de_lfc1_genes"), sub)
  
  write.table(sub, file = paste0("analysis/rna_seq/differential_genes/", as.character(design$group[i]), "_de_lfc1_genes.csv"), 
              sep = ",", row.names = FALSE)
  
  sub <- subset(res, (padj < 0.05))
  
  assign(paste0(as.character(design$group[i]), "_de_lfc0_genes"), sub)
  
  write.table(sub, file = paste0("analysis/rna_seq/differential_genes/", as.character(design$group[i]), "_de_lfc0_genes.csv"), 
              sep = ",", row.names = FALSE)
  
  assign(paste0(as.character(design$group[i]), "_de_all"), res)
}

save(list = c(sapply(design$group, function(x) paste0(x, c("_de_all", "_de_lfc0_genes", "_de_lfc1_genes")))), file = "R_Data/Subtype_DE_Genes.RData")
