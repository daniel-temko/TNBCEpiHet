library(Seurat)

setwd("../")

obj.sc <- load(file = "R_Data/SC_BC_DR.RData")

sample3 <- sc$sample3

save(sample3, file = "R_Data/SC_Sample3.RData")
