library(UpSetR)
setwd("../")

load("R_Data/PRRX1_Targets.RData")
load("R_Data/Gene_Sigs.RData")

temp <- list(bas.up = tnbc_gsl$bas$up,
             mes.up = tnbc_gsl$mes$up,
             lum.up = tnbc_gsl$lum$up,
             mes.rna.tgt.up = prrx1_gsl$mes_rna_targets$up,
             mes.rna.tgt.dn = prrx1_gsl$mes_rna_targets$dn)

pdf("analysis/prrx1_chip_seq/targets/mes_rna_target_tnbc_type_overlap.pdf", width = 11, height = 9)
upset(fromList(temp), 
      order.by = "freq", 
      point.size = 3.5, 
      line.size = 2, 
      mainbar.y.label = "Gene Overlaps", 
      sets.x.label = "Number of Genes", 
      #text.scale = c(1.3, 1.3, 1, 1, 2, 1.5),
      text.scale = c(2.6, 2.6, 2, 2, 4, 3))
dev.off()

temp <- list(bas.up = tnbc_gsl$bas$up,
             mes.up = tnbc_gsl$mes$up,
             lum.up = tnbc_gsl$lum$up,
             hs578.rna.tgt.up = prrx1_gsl$hs578_rna_targets$up,
             hs578.rna.tgt.dn = prrx1_gsl$hs578_rna_targets$dn)

pdf("analysis/prrx1_chip_seq/targets/hs578_rna_target_tnbc_type_overlap.pdf", width = 11, height = 9)
upset(fromList(temp), 
      order.by = "freq", 
      point.size = 3.5, 
      line.size = 2, 
      mainbar.y.label = "Gene Overlaps", 
      sets.x.label = "Number of Genes", 
      #text.scale = c(1.3, 1.3, 1, 1, 2, 1.5),
      text.scale = c(2.6, 2.6, 2, 2, 4, 3))
dev.off()

temp <- list(bas.up = tnbc_gsl$bas$up,
             mes.up = tnbc_gsl$mes$up,
             lum.up = tnbc_gsl$lum$up,
             hs578.tgt = prrx1_gml$hs578_targets)

pdf("analysis/prrx1_chip_seq/targets/hs578_target_tnbc_type_overlap.pdf", width = 5, height = 5)
upset(fromList(temp), 
      order.by = "freq", 
      point.size = 3.5, 
      line.size = 2, 
      mainbar.y.label = "Gene Overlaps", 
      sets.x.label = "Number of Genes", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 1.5))
      #text.scale = c(2.6, 2.6, 2, 2, 4, 3))
dev.off()

#--------------------------------------------------------------------------------------------------------------
# Analyse overlap between type genes and mes targets (hs578 targets)

tt <- sapply(tnbc_gsl, function(x) sapply(x, length))

(hs_cc <- sapply(tnbc_gsl, function(x){
   sapply(x, function(y){
      length(intersect(y, prrx1_gml$hs578_targets))
   })
}))

(hs_oo <- sapply(tnbc_gsl, function(x){
   sapply(x, function(y){
      length(intersect(y, prrx1_gml$hs578_targets)) / length(y)
   })
}))


(tt_oo <- sapply(tnbc_gsl, function(x){
   sapply(x, function(y){
      length(intersect(y, prrx1_gml$tt642_targets)) / length(y)
   })
}))

# sig tests
prop.test(x = hs_cc[1,], n = tt[1,])
pairwise.prop.test(x = hs_cc[1,], n = tt[1,])


#--------------------------------------------------------------------------------------------------------------
# Analyse overlap between type genes and mes targets (mes_targets)

tt <- sapply(tnbc_gsl, function(x) sapply(x, length))

(hs_cc <- sapply(tnbc_gsl, function(x){
   sapply(x, function(y){
      length(intersect(y, prrx1_gml$mes_targets))
   })
}))

(hs_oo <- sapply(tnbc_gsl, function(x){
   sapply(x, function(y){
      length(intersect(y, prrx1_gml$mes_targets)) / length(y)
   })
}))


(tt_oo <- sapply(tnbc_gsl, function(x){
   sapply(x, function(y){
      length(intersect(y, prrx1_gml$tt642_targets)) / length(y)
   })
}))

# sig tests
prop.test(x = hs_cc[1,], n = tt[1,])
pairwise.prop.test(x = hs_cc[1,], n = tt[1,])


