setwd("../")

############################################################################################################
## Load libraries
############################################################################################################

source("../utilities/files.R")
library(minfi)
library(DMRcate)
library(limma)
library(RColorBrewer)

############################################################################################################
## Load data
############################################################################################################

obj.mth <- load("R_Data/Methyl_Data.RData")
pal <- brewer.pal(8,"Dark2")

############################################################################################################
## Probe-based analysis
############################################################################################################

cl_type <- factor(methyl_metadata$tnbc.subtype)
design <- model.matrix(~0+cl_type)
colnames(design) <- levels(cl_type)

fit <- lmFit(methyl_all_m, design)

cont_ids <- c("basal", "luminal", "mesenchymal")
cont_mat <- makeContrasts(basal - 0.5*luminal - 0.5*mesenchymal,
                          luminal - 0.5*basal - 0.5*mesenchymal,
                          mesenchymal - 0.5*basal - 0.5*luminal,
                          levels = design)

fit2 <- contrasts.fit(fit, cont_mat)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))

dmps <- topTable(fit2, num = Inf, coef = 2, genelist = methyl_all_ann)

par(mfrow=c(2,2))
sapply(rownames(dmps)[1:4], function(cpg){
  plotCpg(methyl_all_beta, cpg=cpg, pheno=methyl_metadata$tnbc.subtype, ylab = "Beta values")
})

############################################################################################################
## Region-based analysis
############################################################################################################

message("Finding regions")

for(i in 1:length(cont_ids)){
  message(cont_ids[i])
  my_ann <- cpg.annotate(object = methyl_all_m, datatype="array", what="M", analysis.type = "differential",
                         design=design, contrasts=TRUE, cont.matrix = cont_mat, coef = colnames(cont_mat)[i], 
                         arraytype="450K")
  
  dmrs <- dmrcate(my_ann, lambda=1000, C=2)
  results.ranges <- extractRanges(dmrs)
  results.ranges
  res <- as.data.frame(results.ranges)
  
  # write to file after converting to 0-based half-open
  fn <- paste0("analysis/dna_methylation/dmr_de_novo/", cont_ids[i], "_dmr.csv")
  out <- res
  out$start <- out$start + 1
  write.table(out, fn, row.names = FALSE, sep = ",")
  
  # vizualize results
  groups <- pal[1:length(unique(cl_type))]
  names(groups) <- levels(cl_type)
  cols <- groups[as.character(cl_type)]
  
  # plot top DMR
  results.subranges <- results.ranges[which(results.ranges$meandiff<0),]
  fn <- paste0("analysis/dna_methylation/dmr_de_novo/", cont_ids[i], "_example_region.pdf")
  pdf(fn)
  par(mfrow=c(1,1))
  DMR.plot(ranges = results.subranges, dmr = 1, CpGs = methyl_all_beta, phen.col = cols,
           what = "Beta", arraytype = "450K", genome = "hg19")
  dev.off()
  
  assign(paste0("dmr_", cont_ids[i]), res)
}

save(list = paste0("dmr_", cont_ids), file = "R_Data/Methyl_DMR_De_Novo.RData")
