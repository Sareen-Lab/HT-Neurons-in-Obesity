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

### Subset metadata
df.metadata.b <- df.metadata[!row.names(df.metadata) == 'iHT_02iOBS',]
df.metadata.HT <- df.metadata.b[df.metadata.b$Type == 'HT',]
df.metadata.iHT <- df.metadata.b[df.metadata.b$Group == 'iHT',]

### Subset raw data
df.HT.raw <- df.counts.raw[row.names(df.metadata.HT)]
df.iHT.raw <- df.counts.raw[row.names(df.metadata.iHT)]

### Place OBS before CTR to ensure correct fold change direction
df.iHT.raw <- df.iHT.raw[c(8:11,1:7)]

### Convert to integer matrixes
matrix.HT.raw <- round(as.matrix(df.HT.raw[match(row.names(df.metadata.HT),names(df.HT.raw))]))
matrix.iHT.raw <- round(as.matrix(df.iHT.raw[match(row.names(df.metadata.iHT),names(df.iHT.raw))]))
matrix.iHT.raw <- round(as.matrix(df.iHT.raw))

head(matrix.HT.raw)
head(matrix.iHT.raw)

### Perform DESeq2
de.HT <- DESeqDataSetFromMatrix(countData = matrix.HT.raw, colData = df.metadata.HT, design = ~ Source)
de.iHT <- DESeqDataSetFromMatrix(countData = matrix.iHT.raw, colData = df.metadata.iHT, design = ~ Disease)

# Generate results
(de.HT <- DESeq(de.HT))
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

### Subset for low p-value
low.p.HT <- which(df.de.HT$`iHT-v-aHT_padj`<0.00001)
test <- df.de.HT[low.p.HT,]
test.exprs <- df.HT.raw[low.p.HT,]
(nrow(test))

### Compare means
HT.mean <- apply(test.exprs,1,mean)
iHT.mean <- apply(test.exprs[1:11],1,mean)
aHT.mean <- apply(test.exprs[13:17],1,mean)
#data.frame(test$`iHT-v-aHT_baseMean`,HT.mean,iHT.mean,aHT.mean)[1:24,]

### Compare fold changes
data.frame(iHT.mean,aHT.mean,HT.mean,test$`iHT-v-aHT_baseMean`,test$`iHT-v-aHT_log2FC`,log2(iHT.mean/aHT.mean))[1:20,]
### PASS: Sign and magnitude match results when iHT divided by aHT

### iHT sanity check

### Subset for low p-value
low.p.iHT <- which(df.de.iHT$`Obs-v-Ctr_padj`<0.0001)
test <- df.de.iHT[low.p.iHT,]
test.exprs <- df.iHT.raw[low.p.iHT,]
(nrow(test))

### Compare means
iHT.mean <- apply(test.exprs,1,mean)
obs.mean <- apply(test.exprs[1:4],1,mean)
ctr.mean <- apply(test.exprs[5:11],1,mean)
data.frame(obs.mean,ctr.mean,iHT.mean,test$`Obs-v-Ctr_baseMean`,test$`Obs-v-Ctr_log2FC`,log2(obs.mean/ctr.mean))[1:24,]

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

### Generate Differential expression table
df.diff.ex <- df.composite

### Add Differential Expression 1
df.diff.ex <- cbind(df.diff.ex,df.de.HT)

### Add differential Expression 2
df.diff.ex <- cbind(df.diff.ex,df.de.iHT)

### Save normalized
df.normalized <- cbind(df.composite,df.TPM.normalized)

### Save raw
df.raw <- cbind(df.composite,df.counts.raw)

########################################################################
### Output

### Save expression table
setwd("Z:/Data/RNAseq HT neurons and tissue/Andrews_files/")

write.csv(df.diff.ex,"Sareen-HT-2017 data for GEO - Diff-expr.csv", row.names = FALSE)
write.csv(df.normalized,"Sareen-HT-2017 data for GEO - TPM normalized.csv", row.names = FALSE)
write.csv(df.raw,"Sareen-HT-2017 data for GEO - Raw read counts.csv", row.names = FALSE)


########################################################################
### Minor formatting fixes

file.names <- read.csv('Z:/Data/RNAseq HT neurons and tissue/TEMP/file-names-temp.csv', header = FALSE)[,1]

#reformatted.file.names <- data.frame('first' = character(), 'second' = character(), 'third' = character(), 'fourth' = character(), stringsAsFactors = FALSE)
reformatted.file.names <- matrix(1,1,4)

for (sample in 1:20) {
  first.row <- (sample - 1)*4 + 1
  first = as.character(file.names[first.row])
  second = as.character(file.names[first.row +1])
  third = as.character(file.names[first.row + 2])
  fourth = as.character(file.names[first.row + 3])
  new.row <- c(first, second, third, fourth)
  #print(new.row)
  reformatted.file.names <- rbind(reformatted.file.names, new.row)
}

df.reformatted.file.names <- data.frame(reformatted.file.names)[2:21,]

output <- df.reformatted.file.names[c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13),]

write.csv(output,"z:/Data/RNAseq HT neurons and tissue/Andrews_files/TEMP.csv")
