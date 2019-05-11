# script calculates log2FC by computing medians for amp and non amp samples for every gene
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidyverse)
library(data.table)

# read matrix to get FPKM expression values
matrix <- read.delim("GSE49711_SEQC_NB_TUC_G_log2.txt", header = T)
matrix <- matrix %>%
  gather("Sample_title",value,-c("X00gene_id"))

# splitting sample_titles
setDT(matrix)[, c("V1","V2","V3","V4","V5") := tstrsplit(matrix$Sample_title, "_", fixed=TRUE)]
setDT(matrix)[,"new_sample_title" := paste0(matrix$V1,"_",matrix$V2)]
matrix <-  subset(matrix, select = c("X00gene_id","new_sample_title","value"))
colnames(matrix)[colnames(matrix)=="new_sample_title"] <- "Sample_title"




# reading in soft file to get risk and mycn status for sample_titles
family.soft <-  read.delim("GSE49711_family.soft", header = F)
# Sample_titles
sample.title <-  as.data.frame(family.soft[ (family.soft$V1 %like% "!Sample_title" ), ])
colnames(sample.title) <- "V1"
setDT(sample.title)[, c("V1", "sample_title") := tstrsplit(sample.title$V1, "=", fixed=TRUE)]

# high_risk
risk <-  as.data.frame(family.soft[ (family.soft$V1 %like% "!Sample_characteristics_ch1 = high risk"), ])
colnames(risk) <- "V1"
setDT(risk)[, c("V1", "high_risk") := tstrsplit(risk$V1, ":", fixed=TRUE)]

# mycn_status
mycn_status <- as.data.frame(family.soft[ (family.soft$V1 %like% "!Sample_characteristics_ch1 = mycn status"), ])
colnames(mycn_status) <- "V1"
setDT(mycn_status)[, c("V1", "mycn_status") := tstrsplit(mycn_status$V1, ":", fixed=TRUE)]

# combining sample with mycn status and risk
sample.mycn.risk <- as.data.frame(cbind(sample.title$sample_title, mycn_status$mycn_status, risk$high_risk))
colnames(sample.mycn.risk) <-  c("Sample_title","mycn_status","high_risk")
# Sample titles were of type integer, changing to character
sample.mycn.risk$Sample_title <- as.character(sample.mycn.risk$Sample_title)
# removing blank space before the sample_title
sample.mycn.risk$Sample_title <- gsub(" ","",sample.mycn.risk$Sample_title, fixed = T)



# merging matrix with metadata
matrix.mycn <- merge(matrix,sample.mycn.risk, by = "Sample_title")

# calculating median grouping by genes and mycn_status
computed.values <-  matrix.mycn %>%
  filter(., mycn_status != " N/A") %>%
  group_by(X00gene_id,mycn_status) %>%
  summarise(median = median(value)) %>%
  spread(key = "mycn_status", value = "median")

setnames(computed.values, " 0","non-amplified")
setnames(computed.values, " 1","amplified")

setDT(computed.values)[,log2FC:= log2(computed.values$`non-amplified`/computed.values$amplified)]
  
