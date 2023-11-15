setwd("../")

load("R_Data/HMS_Data.RData")

write.table(hms_norm, 
            file = "analysis/hms/summary/hms_norm.csv", 
            sep = ",",
            col.names = NA)

write.table(hms_rep_norm,
            file = "analysis/hms/summary/hms_rep_norm.csv",
            sep = ",",
            col.names = NA)
