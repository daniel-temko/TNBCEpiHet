setwd("../")
library(DESeq2)

# Load data --------------------------------------------------------------------------------

load("R_Data/Subtype_DE_Genes.RData")
imm_sigs <- read.table("immune_gene_sets/combined_immune_sigs.csv", 
                       sep = ",", 
                       header = TRUE, 
                       stringsAsFactors = FALSE)
ps_genes <- read.table("annotation/mammaprint_simona_20200616.csv",
                       sep = ",",
                       header = TRUE,
                       stringsAsFactors = FALSE)
mbs_genes <- read.table("annotation/werb_49_metastasis_sig.txt",
                        header = TRUE,
                        stringsAsFactors = FALSE)
rt_genes <- read.table("annotation/residual_tumor_sig_balko_2012.csv",
                       sep = ",",
                       header = TRUE,
                       stringsAsFactors = FALSE)
imm_n_genes <- read.table("annotation/nelly_immune_signatures_20210317.csv",
                          sep = ",",
                          header = TRUE,
                          stringsAsFactors = FALSE,
                          row.names = 1)
imm_n_meta <- read.table("metadata/Immune_Signature_Names.csv",
                         sep = ",",
                         header = TRUE,
                         stringsAsFactors = FALSE,
                         row.names = 1)

# Create signatures ------------------------------------------------------------------------------------------------

tnbc_gsl <- lapply(list(basal_de_lfc0_genes, luminal_de_lfc0_genes, mesenchymal_de_lfc0_genes), function(x){
  x <- x[order(abs(x$stat), decreasing = TRUE),]
  list(up = x$gene.symbol[which(x$log2FoldChange > 0)],
       dn = x$gene.symbol[which(x$log2FoldChange < 0)])
})
names(tnbc_gsl) <- c("bas", "lum", "mes")

imm_gsl <- lapply(imm_sigs, function(x) {
  list(gene = x[which(x != "")])
})

sv_gsl <- list(ps = list(gene = ps_genes$Gene),
               mbs = list(gene = mbs_genes$Gene),
               rt = list(gene = rt_genes$gene))

ann_list <- lapply(imm_n_genes$ImmPortCategory, function(x) {
  strsplit(x, ", ")[[1]]
})

categories <- unique(unlist(ann_list))

imm_n_gsl <- list()
for(i in 1:length(categories)){
  ids <- sapply(ann_list, function(x) categories[i] %in% x)
  imm_n_gsl[[categories[i]]] <- list(gene = rownames(imm_n_genes)[ids])
}

# Update signature names to remove special characters
new_names <- gsub("( |-)", "_", names(imm_n_gsl))
names(imm_n_gsl) <- new_names

# Process signature metadata
imm_n_map <- imm_n_meta$Short.Name
names(imm_n_map) <- imm_n_meta$Signature

save(tnbc_gsl, imm_gsl, sv_gsl, imm_n_gsl, file = "R_Data/Gene_Sigs.RData")
save(imm_n_map, file = "R_Data/Gene_Sigs_Metadata.RData")
