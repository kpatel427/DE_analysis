
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager",repos = "http://cran.us.r-project.org", dependencies = T)}
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
BiocManager::install("GenomicFeatures")

library(edgeR)
library(GenomicFeatures)
library(ggplot2)
library(reshape)
library(dplyr)
library(tidyverse)
library(DESeq2)

#read in MYCN status
mycn <- read.delim("dataset-MYCN-status.txt", sep = "\t", header = T)
mycn <- mycn %>% column_to_rownames(var = "Sample_title")
setnames(mycn,"mycn_status","Status")
mycn$Status <- gsub(" Amplified","Amplified", mycn$Status)
mycn$Status <- gsub(" Non Amplified","NonAmplified", mycn$Status)
levels(mycn$Status) <-  c("nonAmplified","Amplified")

counts <- read.delim("dataset-count-matrix.txt", sep = "\t", header = TRUE)

##construct DESeq dataset for MYCN status:
dataset <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = mycn,
                                  design = ~ Status)

##remove rows with 0 reads
dataset <- dataset[ rowSums(counts(dataset)) > 0, ]

#amplified vs non-amp written 'non, amp' means +FC = greater in amp subset
dataset$Status <- factor(dataset$Status, levels=c("NonAmplified", "Amplified"))

#differential expression analysis
dataset <- DESeq(dataset)
res <- results(dataset)
resOrdered <- res[order(res$padj),]

#summary for distribution of DE-genes
summary(resOrdered)

#how many adj p-values were < 0.1?
sum(res$padj < 0.1, na.rm=TRUE)

#plot MA
plotMA(res, main="DESeq2", ylim=c(-2,2))

#reorder by  significance
res.1 <- subset(resOrdered, padj < 0.1)

#save table
write.table(as.data.frame(res.1), paste(Sys.Date(), "-DiffExp_MYCN_adjp.10.txt",
                                        sep = ""), row.names = T, quote = F, sep = "\t")

#for list of genes
diffexpgenes <- sort(as.factor(rownames(res.1)))
