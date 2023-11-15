setwd("../")

load(file = "R_Data/Metabric_Annotated.RData")
load(file = "R_Data/TCGA_Annotated.RData")

write.table(mbc_tnbc_metadata, 
            "analysis/correspondence/mbc_tnbc_metadata.csv",
            sep = ",",
            row.names = FALSE)

write.table(tcga_tnbc_metadata,
            "analysis/correspondence/tcga_tnbc_metadata.csv",
            sep = ",",
            row.names = FALSE)