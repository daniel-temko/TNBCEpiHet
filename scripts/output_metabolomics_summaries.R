# Load metabolomics data into R
library(stringr)

setwd("../")

# Load data and merge batches --------------------------------------------------------------------------------------------------------------------------------------------

load("R_Data/Metab_Metadata.RData")

metab1 <- read.table("metabolomics/042617QbjovanovicMets1_54.txt", sep = "\t", 
                     stringsAsFactors = FALSE, 
                     na.strings = "N/A", 
                     row.names = 1, 
                     header = TRUE)[-(1:2),]

orig_names1 <- colnames(metab1)
sample_nums1 <- sapply(colnames(metab1), function(x) strsplit(x, "MET")[[1]][2])
sample_names1 <- paste0("MET", str_pad(sample_nums1, 3, side = "left", pad = "0"))
colnames(metab1) <- sample_names1

metab2 <- read.table("metabolomics/050517QbjovanovicMets55_108.txt", sep = "\t", 
                     stringsAsFactors = FALSE, 
                     na.strings = "N/A", 
                     row.names = 1, 
                     header = TRUE)[-(1:2),]

orig_names2 <- colnames(metab2)
sample_nums2 <- sapply(colnames(metab2), function(x) strsplit(x, "QBJ")[[1]][2])
sample_names2 <- paste0("MET", str_pad(sample_nums2, 3, side = "left", pad = "0"))
colnames(metab2) <- sample_names2

# write name mapping
name_map <- data.frame(name.in.file = c(orig_names1, orig_names2),
                       name.in.metadata = c(sample_names1, sample_names2))

metab_raw <- apply(cbind(metab1, metab2), 2, as.numeric)
rownames(metab_raw) <- rownames(metab1)

# subset samples to match metadata
metab_raw <- metab_raw[,colnames(metab_raw) %in% metab_sample_table$Sample.Name]

# check columns match metadata
all(colnames(metab_raw) == metab_sample_table$Sample.Name)

# Remove all NA dims, batch normalize, and merge replicates --------------------------------------------------------------------------------------------------------------------------------------------

# remove metabolites with no data
keep <- apply(metab_raw, 1, function(x) length(which(is.na(x))) != ncol(metab_raw))
metab_raw <- metab_raw[keep,]

# batch normalisation
batch_mean1 <- apply(metab_raw[,which(metab_sample_table$Batch == 1)], 1, function(x) mean(x, na.rm = TRUE))
batch_mean2 <- apply(metab_raw[,which(metab_sample_table$Batch == 2)], 1, function(x) mean(x, na.rm = TRUE))
overall_mean <- apply(metab_raw, 1, function(x) mean(x, na.rm = TRUE))
Bx <- metab_raw
for(i in 1:ncol(Bx)){
  if(metab_sample_table$Batch[i] == 1){
    Bx[,i] <- Bx[,i] / (batch_mean1 / overall_mean)
  } else {
    Bx[,i] <- Bx[,i] / (batch_mean2 / overall_mean)
  }
}
metab_bc <- Bx

# combine replicates
Cx <- aggregate(x = t(metab_bc), by = list(metab_sample_table$Cell.Line), FUN = function(x) mean(x, na.rm = TRUE))
rownames(Cx) <- Cx$Group.1
Cx$Group.1 <- NULL
metab_comb <- t(Cx)

# remove any metabolites with remaining NA values
keep <- apply(metab_comb, 1, function(x) length(which(is.na(x))) == 0)
metab_comb <- metab_comb[keep,]

# rename outputs
all(colnames(metab_bc) == metab_sample_table$Sample.Name)
colnames(metab_bc) <- metab_sample_table$Sample.Name2

write.table(name_map, 
            file = "analysis/metabolomics/summary/sample_name_map.csv",
            row.names = FALSE,
            sep = ",")

# write file including all remaining metabolites used in downstream analysis
metab_replicates <- metab_bc[rownames(metab_comb),]
write.table(metab_replicates, 
            file = "analysis/metabolomics/summary/filtered_metabolites_replicates.csv",
            sep = ",",
            col.names = NA)

# write file excluding na metabolites
metab_no_na <- na.omit(metab_replicates)
write.table(metab_no_na, 
            file = "analysis/metabolomics/summary/filtered_metabolites_replicates_na_rm.csv",
            sep = ",",
            col.names = NA)