### Processing RNA-seq expression data for submission to GEO -- Andrew R Gross -- 2017/02/16
### I have a lot of different data in different formats.  I need to collect up the raw counts, the normalized reads, and the diff expression results into one table

########################################################################
### Header

library(DESeq2)
library(biomaRt)

#ensembl = useMart(host="www.ensembl.org")
ensembl = useMart(host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
listMarts(host="www.ensembl.org")
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL')
listDatasets(ensembl)
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
#filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

########################################################################
### Functions

convertIDs <- function(dataframe) {
  ensemblIDs <- c()
  for (rowName in row.names(dataframe)) {
    ensemblID <- strsplit(rowName,"\\.")[[1]][1]
    ensemblIDs <- c(ensemblIDs,ensemblID)
  }
  row.names(dataframe)<-ensemblIDs
  return(dataframe)
}

addDescription <- function(dataframe) {
  descriptions <- getBM(attributes=c('ensembl_gene_id','description'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  #print(descriptions)
  descriptions <- descriptions[match(row.names(dataframe),descriptions[,1]),]
  Descrip <- c()
  for (rowNumber in 1:length(descriptions[,1])) {
    fullDescr <- descriptions[rowNumber,][,2]
    shortDescr <- strsplit(fullDescr," \\[Source")[[1]][1]
    Descrip <- c(Descrip, shortDescr)
  }
  dataframe[length(dataframe)+1] <- Descrip
  names(dataframe)[ncol(dataframe)] <- "Description"
  return(dataframe)
}

addGene <- function(dataframe) {
  genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  genes <- genes[match(row.names(dataframe),genes[,1]),]
  Gene <- c()
  for (rowNumber in 1:length(genes[,1])) {
    newGene <- genes[rowNumber,][,2]
    Gene <- c(Gene, newGene)
  }
  dataframe[length(dataframe)+1] <- Gene
  names(dataframe)[ncol(dataframe)] <- "Gene"
  return(dataframe)
}

########################################################################
### Input

### Metadata
df.metadata <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/samples.csv", row.names=1)

### Normalized
df.TPM.normalized <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)

### Raw, Non-normalized
df.counts.raw <- read.table("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.all.count", row.names=1, header=TRUE)

### Prepared DE files
#deOBSvCTR1 <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/DE_gene_list before filter/OBS-iHT_v_CTR-iHT_DE.csv",row.names=1)
#deOBSvCTR2 <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/DE_gene_list after_filter/OBS_vs_CTR/OBS-iHT_vs_CTR-iHT_filtered_DE_csv_results_REVISED_TEMP.csv",row.names=1)

########################################################################
### Formatting

### Remove decimal from df.counts.raw
df.counts.raw <- convertIDs(df.counts.raw)

### Reorder rows of df.counts.raw
df.counts.raw <- df.counts.raw[c(1,11,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,12,13)]

### Convert sample names for df.counts.raw
data.frame(names(df.counts.raw), row.names(df.metadata))
names(df.counts.raw) <- row.names(df.metadata)

### Reorder rows of df.TPM.norm
df.TPM.normalized <- df.TPM.normalized[c(1,11,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,12,13)]

### Convert sample names for df.TPM.normalized
data.frame(names(df.TPM.normalized), row.names(df.metadata))
#names(df.TPM.normalized) <- row.names(df.metadata)

########################################################################
### Generate DE data

### Subset raw data
df.HT.raw <- df.counts.raw[grep("HT",df.metadata$Type)]
df.iHT.raw <- df.counts.raw[grep("iHT",df.metadata$Group)]

### Subset metadata
df.metadata.b <- df.metadata[!row.names(df.metadata.HT) == 'iHT_02iOBS',]
df.metadata.HT <- df.metadata[df.metadata.b$Type == 'HT',]
df.metadata.iHT <- df.metadata[df.metadata.b$Group == 'iHT',]

### Convert to integer matrixes
matrix.HT.raw <- round(as.matrix(df.HT.raw[match(row.names(df.metadata.HT),names(df.HT.raw))]))
matrix.iHT.raw <- round(as.matrix(df.iHT.raw[match(row.names(df.metadata.iHT),names(df.iHT.raw))]))

head(matrix.HT.raw)
head(matrix.iHT.raw)

### Perform DESeq2
de.HT <- DESeqDataSetFromMatrix(countData = matrix.HT.raw, colData = df.metadata.HT, design = ~ Source)
de.iHT <- DESeqDataSetFromMatrix(countData = matrix.iHT.raw, colData = df.metadata.iHT, design = ~ Disease)

# Generate results
(de.HT <- DESeq(deHT))
(results.de.HT <- results(de.HT))
df.de.HT <- as.data.frame(results.de.HT)
names(df.de.HT)[2] <- 'log2FC'
names(df.de.HT) <- paste0('iHT-v-aHT_',names(df.de.HT))

(de.iHT <- DESeq(de.iHT))
(results.de.iHT <- results(de.iHT))
df.de.iHT <- as.data.frame(results.de.iHT)
names(df.de.iHT)[2] <- 'log2FC'
names(df.de.iHT) <- paste0('Obs-v-Ctr_',names(df.de.iHT))

########################################################################
### Sanity check

### HT sanity check
test <- 
### Subset for low p-value
low.p.HT <- which(df.de.HT$`iHT-v-aHT_padj`<0.0001)
test <- df.de.HT[low.p.HT,]
test.exprs <- df.HT.raw[low.p.HT,]
(nrow(test))

### Compare stats
HT.mean <- apply(test.exprs,1,mean)
iHT.mean <- apply(test.exprs[1:11],1,mean)
aHT.mean <- apply(test.exprs[13:17],1,mean)
data.frame(test$`iHT-v-aHT_baseMean`,HT.mean,iHT.mean,aHT.mean)

########################################################################
### Formatting before joining

### Rename samples in raw counts table
counts.names <- paste0('raw-', names(df.counts.raw))
names(df.counts.raw) <- counts.names

### Rename samples in normalized table
TPM.names <- paste0('norm-', names(df.TPM.normalized))
names(df.TPM.normalized) <- TPM.names

########################################################################
### Joining dataframes

### Define new dataframe to receive all others
df.composite <- data.frame('Ensembl' = row.names(df.TPM.normalized))
row.names(df.composite) <- row.names(df.counts.raw)

### Add gene names and abreviations
df.composite <- addGene(df.composite)
df.composite <- addDescription(df.composite)
names(df.composite)[3] <- 'Description'

### Add Differential Expression 1
df.composite <- cbind(df.composite,df.de.HT)

### Add differential Expression 2
df.composite <- cbind(df.composite,df.de.iHT)

### Add Normalized values
df.composite <- cbind(df.composite,df.TPM.normalized)

### Add raw values
df.composite <- cbind(df.composite,df.counts.raw)

########################################################################
### Output

### Save expression table
setwd("Z:/Data/RNAseq HT neurons and tissue/Andrews_files/")

write.csv(df.composite,"Full processed data for GEO-HT RNAseq.csv", row.names = FALSE)

########################################################################
### Sanity check

test <- head(df.composite,20)
head(test)


### Subset for low p-value
low.p.HT <- which(df.composite$`iHT-v-aHT_padj`<0.0001)
test <- df.composite[low.p.HT,]
test.b <- test[36:54]

### Sanity check HT (iHT v aHT)
iHT <- apply(test.b[1:12],1,mean)
aHT <- apply(test.b[15:19],1,mean)
HT.mean <- apply(test.b[c(1:11,15:19)],1,mean)
log2iHT.aHT <- log2(iHT/aHT)
data.frame(iHT,aHT,log2iHT.aHT,test$'iHT-v-aHT_log2FC',2^(test$'iHT-v-aHT_log2FC'))

### Sanity check OBS v CTR
CTR <- apply(test.b[1:7],1,mean)
OBS <- apply(test.b[8:11],1,mean)
log2CTR.OBS <- log2(OBS/CTR)
data.frame(CTR,OBS,log2CTR.OBS,test$'Obs-v-Ctr_log2FC')

###############################################################################
###############################################################################

########################################################################
### Subsetting
########################################################################

# Raw
iHTcounts <- counts[grep("iHT",metaData$Group)]
df.iHT.raw <- df.counts.raw[grep("iHT",metaData$Group)]
obsHTcounts <- counts[grep("OBS",metaData$Disease)]
ctrHTcounts <- counts[intersect(grep("iHT",metaData$Group),grep("CTR",metaData$Disease))]

########################################################################
### Filtering
########################################################################

metaData.iHT <- metaData[metaData$Group == 'iHT',]
metaData.iHT <- metaData.iHT[c(1,2,3,5,6,7,8,9,10,11),]
countMatrix <- as.matrix(df.HT.raw[match(row.names(df.metadata.HT),names(df.HT.raw))])

countMatrix <- round(countMatrix)
head(countMatrix)

metaData.selected <- metaData.iHT
#metaData.selected <- metaData.iHT.f

######################################################################################
### Perform DESeq2 calculation
######################################################################################

# Run DESeq
deHT <- DESeqDataSetFromMatrix(countData = countMatrix, colData = metaData.selected, design = ~ Disease)

# Generate results
(deHT <- DESeq(deHT))
(res <- results(deHT))

resOrdered <- res[order(res$padj),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
res05 <- results(deHT, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

#mean(unlist(iHTcounts[2,][metaData.iHT$Disease=='OBS']))
#mean(unlist(iHTcounts[2,][metaData.iHT$Disease=='CTR']))
#mean(unlist(iHTcounts[2,]))

#iHTcounts[1,]
#iHTcounts[1,][metaData.iHT$Disease=='CTR']
#iHTcounts[1,][metaData.iHT$Disease=='OBS']


### Plotting
plotMA(res, main="DESeq2", ylim=c(-2,2))
plotCounts(deHT, gene=which.min(res$padj), intgroup="Disease")


### Format for output
resultsDF <- as.data.frame(resOrdered)

### Fold change cutoffs
logFoldChange <- resultsDF$log2FoldChange
absLFC <- abs(logFoldChange)
lowVals <- absLFC<0.33
lowPvals <- resultsDF$padj>0.1
logFoldChange[lowVals] <- 0
logFoldChange[is.na(logFoldChange)] <- 0
logFoldChange[lowPvals] <- 0
resultsDF$log2FoldChange <- logFoldChange

# Reorder orginal counts
#reorderedHT <- iHTcounts[match(row.names(resultsDF),row.names(iHTcounts)),]
#reorderedOBS <- obsHTdata2[match(row.names(resultsDF),row.names(obsHTdata2)),]
#reorderedCTR <- ctrHTdata2[match(row.names(resultsDF),row.names(ctrHTdata2)),]

#head(reorderedOBS)
#head(reorderedCTR)
head(resultsDF)
#head(reorderedHT)
resultsDF <- addGene(resultsDF)
resultsDF <- addDescription(resultsDF)
names(resultsDF) <- c("baseMean","log2FoldChange", "lfcSE","stat","pvalue","padj","Gene","Descrip")

#outputDataFrame <- data.frame(row.names(resultsDF),reorderedOBS$median,reorderedCTR$median,reorderedOBS$sd,reorderedCTR$sd,resultsDF$log2FoldChange,resultsDF$padj)
outputDataFrame <- data.frame(row.names(resultsDF),resultsDF$Gene,resultsDF$log2FoldChange,resultsDF$padj,resultsDF$Descrip)
names(outputDataFrame) <- c("Ensembl-ID","Gene","Log2FoldChange","p-adj","Description")


write.table(as.data.frame(outputDataFrame[1:8000,]), file="z://Data/RNAseq HT neurons and tissue/2nd rerun/DE_genes_ANDREW/Obs_iHT_vs_Ctr_iHT--8000.txt",sep="\t",row.names = FALSE)



