setwd("../")

##########################################################################################
## Load library
##########################################################################################

library(ggplot2)
library(cowplot)
library(ggrepel)
source("../utilities/utils.R")

##########################################################################################
## Local functions
##########################################################################################

plot_expression_comparison <- function(f.pp, repel_force = 50, ymax = 6.7, xmax = 6.7){
  keep <- f.pp$sum > 0 | f.pp$hcc > 0
  pp_s <- f.pp[keep,]
  res <- wilcox.test(pp_s$sum - pp_s$hcc)
  
  p <- ggplot(f.pp, aes(x = hcc, y = sum)) + geom_point() + 
    geom_abline(slope = 1, intercept = 0) + 
    geom_text_repel(aes(label = tf), force = repel_force) +
    ylab("SUM185 Expression") +
    xlab("HCC3153 Expression") +
    ylim(c(0, ymax)) +
    xlim(c(0, xmax)) +
    my_theme() +
    geom_text(x = max(f.pp$hcc) * 0.2, y = max(f.pp$sum) * 0.8, label = paste0("P=", round(res$p.value, 3)))
  return(p)
}

##########################################################################################
## Load data
##########################################################################################

obj.rna <- load("R_Data/RNA_Data.RData")

hoc_ann <- read.table("../ref/homer/motif_databases/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv",
                      sep = "\t", header = TRUE)

hcc_st_up_comm_homer <- read.table("analysis/chip_seq_prrx1_oe/homer/hocomoco_filt_peaks_in_superenhancers_hcc_st_up_comm/knownResults.txt",
                                   sep = "\t", comment.char = "", header = TRUE)
hcc_st_dn_comm_homer <- read.table("analysis/chip_seq_prrx1_oe/homer/hocomoco_filt_peaks_in_superenhancers_hcc_st_dn_comm/knownResults.txt",
                                   sep = "\t", comment.char = "", header = TRUE)
hcc_lt_up_comm_homer <- read.table("analysis/chip_seq_prrx1_oe/homer/hocomoco_filt_peaks_in_superenhancers_hcc_lt_up_comm/knownResults.txt",
                                   sep = "\t", comment.char = "", header = TRUE)
hcc_lt_dn_comm_homer <- read.table("analysis/chip_seq_prrx1_oe/homer/hocomoco_filt_peaks_in_superenhancers_hcc_lt_dn_comm/knownResults.txt",
                                   sep = "\t", comment.char = "", header = TRUE)

hcc_st_up_uniq_homer <- read.table("analysis/chip_seq_prrx1_oe/homer/hocomoco_filt_peaks_in_superenhancers_hcc_st_up_uniq/knownResults.txt",
                                   sep = "\t", comment.char = "", header = TRUE)
hcc_st_dn_uniq_homer <- read.table("analysis/chip_seq_prrx1_oe/homer/hocomoco_filt_peaks_in_superenhancers_hcc_st_dn_uniq/knownResults.txt",
                                   sep = "\t", comment.char = "", header = TRUE)
hcc_lt_up_uniq_homer <- read.table("analysis/chip_seq_prrx1_oe/homer/hocomoco_filt_peaks_in_superenhancers_hcc_lt_up_uniq/knownResults.txt",
                                   sep = "\t", comment.char = "", header = TRUE)
hcc_lt_dn_uniq_homer <- read.table("analysis/chip_seq_prrx1_oe/homer/hocomoco_filt_peaks_in_superenhancers_hcc_lt_dn_uniq/knownResults.txt",
                                   sep = "\t", comment.char = "", header = TRUE)

##########################################################################################
## Transform data
##########################################################################################

# Create long format results data
dd <- rbind(data.frame(group = rep("st_up_comm", nrow(hcc_st_up_comm_homer)), motif = hcc_st_up_comm_homer$Motif.Name, rank = 1:nrow(hcc_st_up_comm_homer), logp = hcc_st_up_comm_homer$Log.P.value, q = hcc_st_up_comm_homer$q.value..Benjamini.),
            data.frame(group = rep("st_dn_comm", nrow(hcc_st_dn_comm_homer)), motif = hcc_st_dn_comm_homer$Motif.Name, rank = 1:nrow(hcc_st_dn_comm_homer), logp = hcc_st_dn_comm_homer$Log.P.value, q = hcc_st_dn_comm_homer$q.value..Benjamini.),
            data.frame(group = rep("lt_up_comm", nrow(hcc_lt_up_comm_homer)), motif = hcc_lt_up_comm_homer$Motif.Name, rank = 1:nrow(hcc_lt_up_comm_homer), logp = hcc_lt_up_comm_homer$Log.P.value, q = hcc_lt_up_comm_homer$q.value..Benjamini.),
            data.frame(group = rep("lt_dn_comm", nrow(hcc_lt_dn_comm_homer)), motif = hcc_lt_dn_comm_homer$Motif.Name, rank = 1:nrow(hcc_lt_dn_comm_homer), logp = hcc_lt_dn_comm_homer$Log.P.value, q = hcc_lt_dn_comm_homer$q.value..Benjamini.),
            data.frame(group = rep("st_up_uniq", nrow(hcc_st_up_uniq_homer)), motif = hcc_st_up_uniq_homer$Motif.Name, rank = 1:nrow(hcc_st_up_uniq_homer), logp = hcc_st_up_uniq_homer$Log.P.value, q = hcc_st_up_uniq_homer$q.value..Benjamini.),
            data.frame(group = rep("st_dn_uniq", nrow(hcc_st_dn_uniq_homer)), motif = hcc_st_dn_uniq_homer$Motif.Name, rank = 1:nrow(hcc_st_dn_uniq_homer), logp = hcc_st_dn_uniq_homer$Log.P.value, q = hcc_st_dn_uniq_homer$q.value..Benjamini.),
            data.frame(group = rep("lt_up_uniq", nrow(hcc_lt_up_uniq_homer)), motif = hcc_lt_up_uniq_homer$Motif.Name, rank = 1:nrow(hcc_lt_up_uniq_homer), logp = hcc_lt_up_uniq_homer$Log.P.value, q = hcc_lt_up_uniq_homer$q.value..Benjamini.),
            data.frame(group = rep("lt_dn_uniq", nrow(hcc_lt_dn_uniq_homer)), motif = hcc_lt_dn_uniq_homer$Motif.Name, rank = 1:nrow(hcc_lt_dn_uniq_homer), logp = hcc_lt_dn_uniq_homer$Log.P.value, q = hcc_lt_dn_uniq_homer$q.value..Benjamini.))

stopifnot(all(dd$motif %in% hoc_ann$Model))

dd$tf <- hoc_ann$Transcription.factor[match(dd$motif, hoc_ann$Model)]

stopifnot(all(dd$tf %in% rownames(rna_rpkm_norm_uf)))

dd$sum <- rna_rpkm_norm_uf[dd$tf, "SUM185"]
dd$hcc <- rna_rpkm_norm_uf[dd$tf, "HCC3153"]

gg <- subset(dd, q < 0.05 & rank <= 10)

table(gg$group)

pdf("analysis/chip_seq_prrx1_oe/putative_cofactors/st_dn_comm_expression.pdf")
pp <- subset(gg, group == "st_dn_comm")
plot_expression_comparison(pp)
dev.off()

pdf("analysis/chip_seq_prrx1_oe/putative_cofactors/lt_up_comm_expression.pdf")
pp <- subset(gg, group == "lt_up_comm")
plot_expression_comparison(pp)
dev.off()

pdf("analysis/chip_seq_prrx1_oe/putative_cofactors/lt_dn_comm_expression.pdf")
pp <- subset(gg, group == "lt_dn_comm")
plot_expression_comparison(pp)
dev.off()

# Combined lt up and down plot
pdf("analysis/chip_seq_prrx1_oe/putative_cofactors/lt_up_and_dn_comm_expression.pdf")
pp1 <- subset(gg, group == "lt_up_comm")
p1 <- plot_expression_comparison(pp1, repel_force = 50)
pp2 <- subset(gg, group == "lt_dn_comm")
p2 <- plot_expression_comparison(pp2, repel_force = 50)
plot_grid(p1, p2, align = "h", ncol = 2)
dev.off()
