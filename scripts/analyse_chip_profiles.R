#####
# create H3K27ac genome track plots at type-specific SE's
#####

setwd("../")

################################################################################
# load library
################################################################################

library(Gviz)
library(rtracklayer)
library(colorspace)
source("../utilities/utils.R")

################################################################################
# load data
################################################################################

obj.chip <- load("R_Data/ChIP_Data.RData")
obj.meta <- load("R_Data/ChIP_Metadata.RData")
obj.heat <- load("R_Data/Heatmap_Metadata.RData")

types <- c("basal", "luminal", "mesenchymal")
chip_treat_samp <- subset(chip_sample_table, Pulldown=="H3K27ac")

n_samp <- 3

rids <- c("chr1:170630495-170684023", "chr14:38013196-38096409", "chr5:1881970-1890803")

# get pal
ss <- chip_treat_samp[order(chip_treat_samp$Cell.Line),]

top_list <- lapply(rids, function(x){
    top_lines <- lapply(types, function(y){
        xx <- subset(ss, tnbc.subtype == y)
        stopifnot((x %in% rownames(chip_se_all_rpkm_norm)) & 
                  all(xx$Cell.Line %in% colnames(chip_se_all_rpkm_norm)))
        names(sort(chip_se_all_rpkm_norm[x, xx$Cell.Line], decreasing = TRUE))[1:n_samp]
    })
    names(top_lines) <- types
    top_lines
})

# subset to only lines to be plotted
ss.pal <- subset(ss, Cell.Line %in% unlist(top_list))

pal <- unlist(unname(sapply(types, function(x){
    xx <- subset(ss.pal, tnbc.subtype == x)
    base_col <- line_colours$type[x]
    top_col <- darken(base_col, 0.5)
    col_vals <- colorRampPalette(c(base_col, top_col))(nrow(xx) + 1)[-1]
    names(col_vals) <- xx$Cell.Line
    col_vals
})))

################################################################################
# local function
################################################################################

profile_plot <- function(rid, stacking, f.pal = pal, f.ss = ss, dd = chip_se_all_rpkm_norm, n_sample=n_samp, y_max = 5, gen="hg19", tt=types, sp="chip_seq_cell_lines_cobra", sf="_ds0_treat_pileup.bw"){
    gr <- granges_from_region_ids(rid)
    chr <- as.character(seqnames(gr))
    lm <- round(extendrange(r=c(start(gr), end(gr)), f=0.1))
    
    # get top lines
    top_lines <- lapply(tt, function(x){
        xx <- subset(f.ss, tnbc.subtype == x)
        stopifnot((rid %in% rownames(dd)) & 
                  all(xx$Cell.Line %in% colnames(dd)))
        names(sort(dd[rid, xx$Cell.Line], decreasing = TRUE))[1:n_sample]
    })
    names(top_lines) <- tt
    
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome = gen, chromosome = chr)
    
    tracks <- lapply(tt, function(x){
        xx <- subset(f.ss, Cell.Line %in% top_lines[[x]])
        tl <- lapply(1:nrow(xx), function(y){
            fp <- file.path(sp, tolower(xx$Sample.Name[y]))
            fn <- paste0(xx$prefix[y], sf)
            ff <- file.path(fp, fn)
            dt <- DataTrack(range = ff, genome = gen, type = "l", chromosome = chr, 
                            col = f.pal[sort(top_lines[[x]])], groups = factor(xx$Cell.Line[y],
                            levels = sort(top_lines[[x]])), legend = TRUE)
        })
        ot <- OverlayTrack(tl)
    })
    refGenes <- UcscTrack(genome = gen, chromosome = chr, track = "NCBI RefSeq", table = "ncbiRefSeq", 
                          from = lm[1], to = lm[2], trackType = "GeneRegionTrack", rstarts = "exonStarts",
                          rends = "exonEnds", gene = "name", symbol = "name2", transcript = "name",
                          strand = "strand", fill = "#8282d2", name = "RefSeq genes", stacking = stacking,
                          transcriptAnnotation = "symbol")
    ht <- HighlightTrack(trackList = c(tracks, refGenes), start = start(gr), end = end(gr), chromosome = chr)
    plotTracks(c(gtrack, itrack, ht), from = lm[1], to = lm[2], ylim = c(0, y_max))
    
    # subset refGenes to region
    subset(refGenes@range, start >= lm[1] & end <= lm[2])
}

################################################################################
# plot selected regions
################################################################################

pdf("analysis/chip_seq/gene_tracks/prrx1_squish_plot.pdf")
(gg <- profile_plot(rids[1], "squish", y_max = 8))
dev.off()

pdf("analysis/chip_seq/gene_tracks/prrx1_dense_plot.pdf")
(gg <- profile_plot(rids[1], "dense", y_max = 8))
dev.off()

write.table(gg, "analysis/chip_seq/gene_tracks/prrx1_genes.csv", quote=FALSE, sep=",")

pdf("analysis/chip_seq/gene_tracks/foxa1_squish_plot.pdf")
(gg <- profile_plot(rids[2], "squish", y_max = 8))
dev.off()

pdf("analysis/chip_seq/gene_tracks/foxa1_dense_plot.pdf")
(gg <- profile_plot(rids[2], "dense", y_max = 8))
dev.off()

write.table(gg, "analysis/chip_seq/gene_tracks/foxa1_genes.csv", quote=FALSE, sep=",")

pdf("analysis/chip_seq/gene_tracks/irx4_squish_plot.pdf")
(gg <- profile_plot(rids[3], "squish", y_max = 8))
dev.off()

pdf("analysis/chip_seq/gene_tracks/irx4_dense_plot.pdf")
(gg <- profile_plot(rids[3], "dense", y_max = 8))
dev.off()

write.table(gg, "analysis/chip_seq/gene_tracks/irx4_genes.csv", quote=FALSE, sep=",")
