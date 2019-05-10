library(edgeR)
library(tidyverse)
library(dplyr)
library(data.table)

#read in MYCN status
mycn <- read.delim("dataset-MYCN-status.txt", sep = "\t", header = T)
mycn <- mycn %>% column_to_rownames(var = "Sample_title")
setnames(mycn,"mycn_status","Status")
mycn$Status <- gsub(" Amplified","Amplified", mycn$Status)
mycn$Status <- gsub(" Non Amplified","NonAmplified", mycn$Status)
mycn$Status <- gsub(" NonAmplified","NonAmplified", mycn$Status)
mycn$Status <- as.factor(mycn$Status)


# reading in the count data
counts <- read.delim("data/dataset-count-matrix.txt", sep = "\t")
new.counts <- counts[c(1:60778),c(1,10:507)]

# calculating means for every row grouping by the gene; and selecting rows with max mean
new.counts <- setDT(new.counts)[, .SD[which.max(rowMeans(.SD))], by=`X.Gene`]
new.counts <-  new.counts %>%
  gather(.,"Sample_title",value, -c("X.Gene"))

# back transforming from log base 2
new.counts$value <- 2^(new.counts$value)
new.counts$value <-  as.integer(new.counts$value)

new.counts <- new.counts %>%
  spread(.,Sample_title,value) 

# make gene column into rownames
#rownames(new.counts) <- counts[,1]
new.counts <-  new.counts %>% column_to_rownames(var= "X.Gene")

# In order to keep only those columns in new.counts whose status info is present in new.mycn
new.counts <- new.counts[colnames(new.counts) %in% rownames(mycn)]

# creating a DGE data object
d <- DGEList(counts=new.counts,group=factor(mycn$Status))
d$samples

# we must have at least 100 counts per million (calculated with cpm() in R) 
# on any particular gene that we want to keep
dim(d) # [1] 60586   493
head(cpm(d))

# to calculate counts per million per sample
apply(d$counts, 2, sum) # total gene counts per sample

# keeping genes having count more than 100
keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d) #[1] 6157  493

# resetting library sizes after filtering
d$samples$lib.size <- colSums(d$counts)
d$samples


# normalization
d <- calcNormFactors(d)
d$counts


# dispersion
d1 <- estimateCommonDisp(d, verbose=T)
# Disp = 0.587 , BCV = 0.7662 
d1 <- estimateTagwiseDisp(d1)

# plotting dispersion
plotBCV(d1)

# DE analysis
et12 <- exactTest(d1, pair=c("Amplified","NonAmplified")) 
#topTags(et12, n=6157)

# saving it in a dataframe
tab_DE_genes <- as.data.frame(et12$table)

# making rownames into column
setDT(tab_DE_genes, keep.rownames = TRUE)[]
setnames(tab_DE_genes, "rn","Gene")

# The total number of differentially expressed genes at FDR< 0:05 is:
de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
summary(de1)
tab_genes_values <- as.data.frame(de1)

# merging into one dataframe
final.tab.DE.genes <- cbind(tab_DE_genes,tab_genes_values)

# writing it out to a file
write.table(final.tab.DE.genes, file = paste0(Sys.Date(),"-edgeR-DE-genes-GSE47911-mycnstatus.txt"),
            sep = "\t",quote = F, row.names = F, col.names = T)

# differentially expressed tags from the naive method in d1
de1tags12 <- rownames(d1)[as.logical(de1)] 
plotSmear(et12, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")
