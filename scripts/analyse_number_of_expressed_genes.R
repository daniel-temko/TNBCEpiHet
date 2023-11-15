library(ggplot2)
library(FSA)

setwd("../")

obj.expr <- load("R_Data/RNA_Data.RData")
obj.heat <- load("R_Data/Heatmap_Metadata.RData")

cc <- apply(rna_rpkm, 2, function(x) length(which(x > 1)))
boxplot(cc ~ rna_metadata$tnbc.subtype)
res_kw <- kruskal.test(cc ~ rna_metadata$tnbc.subtype)
res_dn <- dunnTest(cc ~ factor(rna_metadata$tnbc.subtype), method = "holm")

pp <- data.frame(tnbc.type = rna_metadata$tnbc.subtype, ng = cc)
pp$tnbc.type <- factor(pp$tnbc.type, levels = c("luminal", "basal", "mesenchymal"))

pa <- data.frame(tnbc.type = levels(pp$tnbc.type),
                 n = paste0("N=", as.numeric(table(pp$tnbc.type))),
                 x = 1:length(levels(pp$tnbc.type)),
                 y = 13000)

pdf("analysis/rna_seq/diversity/expressed_genes.pdf")
p <- ggplot(pp, aes(x = tnbc.type, y = ng, fill = tnbc.type)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, show.legend = FALSE) +
  scale_fill_manual(values = line_colours$type) +
  ylab("Expressed Genes") +
  theme(text = element_text(size = 18)) +
  geom_text(size = 6,
            data = pa,
            show.legend = FALSE,
            aes(x = x, y = y, label = n))
print(p)
dev.off()

p_val <- res_kw$p.value
write.table(data.frame(p_val = p_val),
            "analysis/rna_seq/diversity/expressed_genes_kruskal_test_p_val.csv",
            row.names = FALSE,
            sep = ",")

write.table(res_dn[[2]],
            "analysis/rna_seq/diversity/expressed_genes_pairwise_comparisons.csv",
            sep = ",",
            col.names = NA)