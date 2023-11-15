load("R_Data/DS_Data.RData")

ds_data <- 1 - ds_auc
xx <- as.matrix(sweep(ds_data, 1, rowMeans(ds_data)))
which(xx == sort(as.matrix(xx), decreasing = T)[1], arr.ind = T)
yy <- as.matrix(sweep(ds_auc, 1, rowMeans(ds_auc)))

which(yy == sort(as.matrix(yy))[1], arr.ind = T)
which(yy == sort(as.matrix(yy))[2], arr.ind = T)
which(yy == sort(as.matrix(yy))[3], arr.ind = T)

vv <- apply(ds_auc, 1, var)
sort(vv, decreasing = T)[1:3]
