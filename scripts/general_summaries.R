setwd("../")

library(reshape)

obj.ds <- load("R_Data/DS_Data.RData")
obj.hms <- load("R_Data/HMS_Data.RData")

## Drug screen results prioritization
xx <- sweep(ds_auc, 1, rowMeans(ds_auc))
dd <- melt(cbind(xx, drug = rownames(xx)), id.vars = "drug")
oo <- melt(cbind(ds_auc, drug = rownames(ds_auc)), id.vars = "drug")
ss <- dd[order(dd$value), ]
ss$orig <- oo[order(dd$value),]$value
ss$av <- rowMeans(ds_auc)[ss$drug]

head(ss)

## HMS expression of cluster-driving marks
dd <- sweep(hms_norm, 1, rowMeans(hms_norm))

ii1 <- grep("H4\\(20to23\\)K20me3", rownames(dd))
rownames(dd)[ii1]
xx1 <- dd[ii1,]
(mm1 <- aggregate(as.numeric(xx1), by = list(hms_metadata$tnbc.subtype), mean))

ii2 <- grep("H3K27ac", rownames(dd))
rownames(dd)[ii2]
xx2 <- colMeans(dd[ii2,])
(mm2 <- aggregate(xx2, by = list(hms_metadata$tnbc.subtype), mean))
