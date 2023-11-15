library(ggplot2)

setwd("../")

# Load data -------------------------------------------------------------------------------------------

load("R_Data/RNA_PRRX1_Data.RData")

# PRRX1 knockdown -------------------------------------------------------------------------------------------

Data <- data.frame(prrx1.rpkm = as.numeric(rna_prrx1_rpkm["PRRX1",]), 
                   prrx1.log.rpkm = as.numeric(rna_prrx1_rpkm_norm["PRRX1",]))
Data$condition <- factor(rna_prrx1_metadata$condition, levels = c("no_DOX", "DOX"))
Data$cell.line <- with(rna_prrx1_metadata, paste(parental.cell.line, shrna))

pdf(file = "analysis/rna_seq_prrx1_kd/comparison/prrx1_rpkm.pdf", width = 10, height = 8)
  p <- ggplot(Data, aes(x = cell.line, y = prrx1.rpkm, fill = condition)) +
          geom_bar(stat = "identity", position=position_dodge()) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
dev.off()

pdf(file = "analysis/rna_seq_prrx1_kd/comparison/prrx1_log_rpkm.pdf", width = 10, height = 8)
  p <- ggplot(Data, aes(x = cell.line, y = prrx1.log.rpkm, fill = condition)) +
          geom_bar(stat = "identity", position=position_dodge()) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
dev.off()


# DDR2 -------------------------------------------------------------------------------------------

Data <- data.frame(ddr2.rpkm = as.numeric(rna_prrx1_rpkm["DDR2",]),
                   ddr2.log.rpkm = as.numeric(rna_prrx1_rpkm_norm["DDR2",]))
Data$condition <- factor(rna_prrx1_metadata$condition, levels = c("no_DOX", "DOX"))
Data$cell.line <- with(rna_prrx1_metadata, paste(parental.cell.line, shrna))

pdf(file = "analysis/rna_seq_prrx1_kd/comparison/ddr2_rpkm.pdf", width = 10, height = 8)
  p <- ggplot(Data, aes(x = cell.line, y = ddr2.rpkm, fill = condition)) +
          geom_bar(stat = "identity", position=position_dodge()) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
dev.off()

pdf(file = "analysis/rna_seq_prrx1_kd/comparison/ddr2_log_rpkm.pdf", width = 10, height = 8)
  p <- ggplot(Data, aes(x = cell.line, y = ddr2.log.rpkm, fill = condition)) +
          geom_bar(stat = "identity", position=position_dodge()) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
dev.off()


# DDR1 -------------------------------------------------------------------------------------------

Data <- data.frame(ddr1.rpkm = as.numeric(rna_prrx1_rpkm["DDR1",]),
                   ddr1.log.rpkm = as.numeric(rna_prrx1_rpkm_norm["DDR1",]))
Data$condition <- factor(rna_prrx1_metadata$condition, levels = c("no_DOX", "DOX"))
Data$cell.line <- with(rna_prrx1_metadata, paste(parental.cell.line, shrna))

pdf(file = "analysis/rna_seq_prrx1_kd/comparison/ddr1_rpkm.pdf", width = 10, height = 8)
  p <- ggplot(Data, aes(x = cell.line, y = ddr1.rpkm, fill = condition)) +
          geom_bar(stat = "identity", position=position_dodge()) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
dev.off()

pdf(file = "analysis/rna_seq_prrx1_kd/comparison/ddr1_log_rpkm.pdf", width = 10, height = 8)
  p <- ggplot(Data, aes(x = cell.line, y = ddr1.log.rpkm, fill = condition)) +
          geom_bar(stat = "identity", position=position_dodge()) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
dev.off()

# IGFBP7 -------------------------------------------------------------------------------------------

Data <- data.frame(igfbp7.rpkm = as.numeric(rna_prrx1_rpkm["IGFBP7",]),
                   igfbp7.log.rpkm = as.numeric(rna_prrx1_rpkm_norm["IGFBP7",]))
Data$condition <- factor(rna_prrx1_metadata$condition, levels = c("no_DOX", "DOX"))
Data$cell.line <- with(rna_prrx1_metadata, paste(parental.cell.line, shrna))

pdf(file = "analysis/rna_seq_prrx1_kd/comparison/igfbp7_rpkm.pdf", width = 10, height = 8)
p <- ggplot(Data, aes(x = cell.line, y = igfbp7.rpkm, fill = condition)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p)
dev.off()

pdf(file = "analysis/rna_seq_prrx1_kd/comparison/igfbp7_log_rpkm.pdf", width = 10, height = 8)
p <- ggplot(Data, aes(x = cell.line, y = igfbp7.log.rpkm, fill = condition)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p)
dev.off()


