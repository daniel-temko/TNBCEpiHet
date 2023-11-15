source("../utilities/diversity.R")
library(ggplot2)

setwd("../")

load("R_Data/RNA_Data.RData")
load("R_Data/Heatmap_Metadata.RData")

rna_sei <- local.rnaseq.shannon(rna_rpkm_uf)

plot_data <- data.frame(tnbc.type = rna_metadata$tnbc.subtype, sei = rna_sei)
plot_data$tnbc.type <- factor(plot_data$tnbc.type, levels = c("luminal", "basal", "mesenchymal"))

plot_ann <- data.frame(tnbc.type = levels(plot_data$tnbc.type),
                       n = paste0("N=", as.numeric(table(plot_data$tnbc.type))),
                       x = 1:length(levels(plot_data$tnbc.type)),
                       y = 0.85)

pdf("analysis/rna_seq/diversity/shannon_equitability.pdf")
  p <- ggplot(plot_data, aes(x = tnbc.type, y = sei, fill = tnbc.type)) + 
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, show.legend = FALSE) + 
    scale_fill_manual(values = line_colours$type) +
    ylab("Shannon's equitability") +
    theme(text = element_text(size = 18)) +
    geom_text(size = 6,
              data = plot_ann,
              show.legend = FALSE,
              aes(x = x, y = y, label = n))
  print(p)
dev.off()

p_val <- kruskal.test(rna_sei, rna_metadata$tnbc.subtype)$p.value
write.table(data.frame(p_val = p_val),
            "analysis/rna_seq/diversity/shann_equitability_kruskal_test_p_val.csv",
            row.names = FALSE,
            sep = ",")