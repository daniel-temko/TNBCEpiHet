setwd("../")

load("R_Data/Metabric_Annotated.RData")

out <- mbc_tnbc_metadata[, c("id2", "tnbc.type"), drop = FALSE]

write.table(out, "analysis/other_analyses/metabric_tnbc_types_for_simona.csv", sep = ",", row.names = FALSE)
