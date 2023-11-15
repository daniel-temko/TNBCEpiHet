source("../utilities/utils.R")
setwd("../")

library(survival)
library(survminer)

tma_data <- read.table("tma/tma_h4k20me3.csv", header = TRUE, sep = ",", na.strings = "ND")
tma_surv <- read.table("tma/tma_survival.csv", header = TRUE, sep = ",")

tma_surv$os.time <- tma_surv$Time.To.Death
tma_surv$os.status <- sapply(tma_surv$Vital.Status..Dead.or.Alive., function(x){
  if(x == "Alive"){
    0
  } else if(x == "Dead"){
    1
  } else {
    NA
  }
})

tma_surv$rfs.time <- tma_surv$Time.To.Recurrence
tma_surv$rfs.status <- sapply(tma_surv$Recurrence, function(x){
  if(x == "No"){
    0
  } else if(x == "Yes"){
    1
  } else{
    NA
  }
})

fit.os <- survfit(Surv(os.time, os.status) ~ H4K20me3.group, data = tma_surv)
pdf("analysis/tma/km_os.pdf", width = 4, height = 4)
ggsurvplot(fit.os, data = tma_surv, pval = TRUE, pval.method = TRUE,
           risk.table = TRUE, ggtheme = my_theme(text_size = 10),
           pval.size = 3, fontsize = 3, tables.height = 0.3)
dev.off()

fit.rfs <- survfit(Surv(rfs.time, rfs.status) ~ H4K20me3.group, data = tma_surv)
pdf("analysis/tma/km_rfs.pdf", width = 4, height = 4)
ggsurvplot(fit.rfs, data = tma_surv, pval = TRUE, pval.method = TRUE,
           risk.table = TRUE, ggtheme = my_theme(text_size = 10),
           pval.size = 3, fontsize = 3, tables.height = 0.3)
dev.off()

res.os <- coxph(Surv(os.time, os.status) ~ H4K20me3.average, data = tma_surv)
res.rfs <- coxph(Surv(rfs.time, rfs.status) ~ H4K20me3.average, data = tma_surv)
