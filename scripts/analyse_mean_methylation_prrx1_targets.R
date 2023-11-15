source("../utilities/utils.R")
source("../utilities/signatures.R")
source("../utilities/plots.R")

mean_m_vals <- function(m_val_data, test_gsl){
  test_gsl_sub <- subset_gsl(test_gsl, rownames(m_val_data))
  mean_m_vals <- as.data.frame(easy_gene_set_score(m_val_data, test_gsl_sub))
}

setwd("../")

# Load data -----------------------------------------------------------------------------------------

load("R_Data/PRRX1_Targets.RData")
load("R_Data/Methyl_Data.RData")
load("R_Data/Gene_SE_Data.RData")
load("R_Data/Heatmap_Metadata.RData")

plot_dir <- "analysis/dna_methylation/mean_methylation/plots"
list_dir <- "analysis/dna_methylation/mean_methylation/lists"

if(!dir.exists(plot_dir)) dir.create(plot_dir)
if(!dir.exists(list_dir)) dir.create(list_dir)

test_gsl <- list(mes_targets = list(gene = prrx1_gml$mes_targets),
                 mes_rna_targets_up = list(gene = prrx1_gsl$mes_rna_targets$up),
                 mes_rna_targets_dn = list(gene = prrx1_gsl$mes_rna_targets$dn))

cut_data_list <- list(se = methyl_gene_se_agg_m,
                      gb = methyl_gb_agg_m,
                      tss = methyl_tss_agg_m)

cut_titles <- data.frame(title = c("Super-enhancers", "Gene Bodies", "TSS"),
                         type = c("superenhancer", "gene_body", "tss"),
                         id = names(cut_data_list))
gs_titles <- data.frame(title = c("Mes Targets", "Mes RNA Targets (PRRX1 High genes)", "Mes RNA Targets (PRRX1 Low genes)"),
                        type = c("mes_targets", "mes_rna_targets_up", "mes_rna_targets_dn"),
                        id = names(test_gsl))

# Analysis -----------------------------------------------------------------------------------------

# Loop over cuts of M-values (with counter i) and gene sets (with counter j)
p_val_list <- list()
for(i in 1:length(cut_titles$id)){
  gs_av_m <- mean_m_vals(cut_data_list[[i]], test_gsl)
  
  # Stat tests
  ## ANOVA
  gs_aov_list <- lapply(gs_av_m, function(x){
    aov(x ~ methyl_metadata$tnbc.subtype)
  })
  gs_p <- sapply(gs_aov_list, function(x) summary(x)[[1]][1, 5])
  p_vals <- data.frame(data_set = cut_titles$id[i], gene_set = names(gs_p), p_val = gs_p)
  p_val_list[[cut_titles$id[i]]] <- p_vals
  
  ## Tukey tests
  gs_tukey_list <- lapply(gs_aov_list, TukeyHSD)
  
  for(j in 1:length(gs_tukey_list)){
    tukey_fname <- paste0(list_dir, "/", cut_titles$type[i], "_", gs_titles$type[j], "_tukey_test.csv")
    write.table(gs_tukey_list[[j]][[1]], tukey_fname, sep = ",", col.names = NA)
  }
}

p_vals <- do.call(rbind, p_val_list)
#p_vals$padj <- p.adjust(p_vals$p_val, metho = "fdr")
write.table(p_vals, paste0(list_dir, "/anova_results.csv"),
            sep = ",", row.names = FALSE)

# Plots
for(i in 1:length(cut_titles$id)){
  gs_av_m <- mean_m_vals(cut_data_list[[i]], test_gsl)
  
  for(j in 1:length(gs_av_m)){
    row_num <- which((p_vals$data_set == cut_titles$id[i]) & (p_vals$gene_set == gs_titles$id[j]))
    #p_adj <- p_vals$padj[row_num]
    p_val <- p_vals$p_val[row_num]
    ylab <- paste0("Average M-value across ", gs_titles$title[j], " in ", cut_titles$title[i])
    
    pdf(paste0(plot_dir, "/", cut_titles$type[i], "_m_vals_in_", gs_titles$type[j], "_by_group.pdf"))
    p <- plot_violin_by_group(gs_av_m[[j]], methyl_metadata$tnbc.subtype, ylab, p_val)  
    print(p)
    dev.off()
  }
}
