
setwd("../")

############################################################################################################
## Load data
############################################################################################################

obj.chip <- load("R_Data/ChIP_Data.RData")
obj.rank <- load("R_Data/ChIP_SE_Ranks.RData")

############################################################################################################
## Local functions
############################################################################################################

rank_ses <- function(dd){
  uu <- lapply(1:ncol(dd), function(x) unique(dd[,x][which(dd[,x]!="")]))
  xx <- unique(unlist(uu))
  nc <- sapply(xx, function(x){
    sum(sapply(uu, function(y) x %in% y))
  })
  rr <- sapply(xx, function(x){
    mean(sapply(uu, function(y) match(x, y)), na.rm = TRUE)
  })
  nc[order(nc, -rr, decreasing = TRUE)]
}

############################################################################################################
## Transform data
############################################################################################################

xx <- se_ranks[,-1]
stopifnot(all(colnames(xx) == rownames(chip_metadata)))

n_lum <- length(which(chip_metadata$tnbc.subtype == "luminal"))
n_bas <- length(which(chip_metadata$tnbc.subtype == "basal"))
n_mes <- length(which(chip_metadata$tnbc.subtype == "mesenchymal"))

keep <- chip_metadata$tnbc.subtype == "luminal"
rr_lum <- rank_ses(xx[,keep])
ff_lum <- rr_lum[(which(rr_lum == sum(keep)))]
pp_lum <- length(ff_lum) / length(rr_lum)

keep <- chip_metadata$tnbc.subtype == "basal"
rr_bas <- rank_ses(xx[,keep])
ff_bas <- rr_bas[(which(rr_bas == sum(keep)))]
pp_bas <- length(ff_bas) / length(rr_bas)

keep <- chip_metadata$tnbc.subtype == "mesenchymal"
rr_mes <- rank_ses(xx[,keep])
ff_mes <- rr_mes[(which(rr_mes == sum(keep)))]
pp_mes <- ff_mes / length(rr_mes)

ff_lum[which(!(names(ff_lum) %in% c(names(ff_bas), names(ff_mes))))]
ff_bas[which(!(names(ff_bas) %in% c(names(ff_lum), names(ff_mes))))]
ff_mes[which(!(names(ff_mes) %in% c(names(ff_bas), names(ff_lum))))]

pdf("analysis/chip_seq/superenhancers/superenhancer_recurrence.pdf")
plot((1:length(rr_lum)) / length(rr_lum), rr_lum / n_lum, col = "blue", ylim = c(0, 1),
     xlab = "Proportion Genes", ylab = "Recurrence")
par(new = TRUE)
plot((1:length(rr_bas)) / length(rr_bas), rr_bas / n_bas, col = "red", ylim = c(0, 1), axes = FALSE,
     xlab = "", ylab = "")
par(new = TRUE)
plot((1:length(rr_mes)) / length(rr_mes), rr_mes / n_mes, col = "green", ylim = c(0, 1), axes = FALSE,
     xlab = "", ylab = "")
dev.off()

