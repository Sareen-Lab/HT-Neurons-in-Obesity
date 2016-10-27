### RNA Seq Tissue comparison by Spearman Correlation, Andrew R Gross, 2016-05-16
### Input: Normalized RNA-seq data; Reference expression data; metadata
### Output: Heatmap displaying pairwise similarity as calculated by Spearman correlation

########################################################################
### Header
########################################################################
library(gplots)
library(RColorBrewer)
library(colorRamps)
library(biomaRt)
library(DESeq2)
library(Hmisc)
listMarts()
ensembl = useMart(host="www.ensembl.org")
ensembl = useMart(host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
listMarts(host="www.ensembl.org")
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL')
listDatasets(ensembl)
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
#filters = listFilters(ensembl)
#attributes = listAttributes(ensembl)

########################################################################
### Functions
########################################################################
addMedSD <- function(dataframe) {   # Add median and standard dev. columns at end of dataframe
  median <- apply(dataframe,1,median)
  sd <- apply(dataframe,1,sd)
  return(data.frame(dataframe,median,sd))  
}
sortByMed <- function(dataframe) {  # Reorder the dataframe by the value in the column 'median'
  order <- order(dataframe$median,decreasing=TRUE)
  return(dataframe[order,])
}
convertIDs <- function(dataframe) { # Remove everything after the decimal in row names
  ensemblIDs <- c()
  for (rowName in row.names(dataframe)) {
    ensemblID <- strsplit(rowName,"\\.")[[1]][1]
    ensemblIDs <- c(ensemblIDs,ensemblID)
  }
  row.names(dataframe)<-ensemblIDs
  return(dataframe)
}
addGene <- function(dataframe) {    # Look up the gene symbol for all the ensembl ids
  genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  genes <- genes[match(row.names(dataframe),genes[,1]),]
  Gene <- c()
  for (rowNumber in 1:length(genes[,1])) { # Reorder the list of gene symbols to match the order of the IDs
    newGene <- genes[rowNumber,][,2]
    Gene <- c(Gene, newGene)
  }
  dataframe[length(dataframe)+1] <- Gene
  names(dataframe)[ncol(dataframe)] <- "Gene"
  return(dataframe)
}
########################################################################
### Import Data
########################################################################
### Normalized
TPMdata <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)
metaData <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/samples.csv", row.names=1)

### Import references
references <- read.table("c://Users/grossar/Bioinform/DATA/rna_seq/Reference_transcriptomes/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct",sep="\t",header=TRUE,row.names=1)

### Reference names, formatted
referenceNames <- c("Adipose, subcutaneous","Adipose, omentum","Adrenal gland","Aorta","Coronary artery","Tibial artery","Bladder","Amygdala","Anteriror cingulate nucleous","Caudate nucleous",
                    "Cerebellar hemisphere","Cerebellum","Cortex","Frontal cortex BA9","Hippocampus","Hypothalamus","Nucleus accumbens","Putamen",
                    "Spinal cord","Substantia nigra","Mammary","Lymphocyte","Fibroblast","Ectocervix","Endocervix","Colon, sigmoid","Colon, transverse","Gastroesophageal junction","Esophagus, mucosa",
                    "Esophagus, muscularis","Fallopian tube","Heart, Atrial","Heart, left ventricle","Kidney","Liver","Lung","Salvitory gland","Skeletal muscle","Tibial nerve","Ovary","Pancreas",
                    "Pituitary","Prostate","Skin, sun-hidden","Skin, sun-exposed","Small intestine","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole blood")

########################################################################
### Formating
########################################################################

### Replace row and column names
sampleNames <- c("iHT_03iCTR","iHT_90iOBS","iHT_77iOBS","iHT_02iOBS","iMN_87iCTR","iMN_201iCTR","aHT_1662_S6","aHT_1838_S7","aHT_1843_S8","aHT_2266_S9","iHT_02iCTR_S13","aHT_2884_S11","21-Reference","iHT_87iCTR","iHT_201iCTR","iHT_25iCTR_S16","iHT_688iCTR_","iHT_80iCTR","iHT_74iOBS","iHT_03iOBS")
names(TPMdata) <- sampleNames
TPMdata <- convertIDs(TPMdata)

### Remove 02iOBS and both iMN
TPMdata[c("iHT_02iOBS","iMN_87iCTR","iMN_201iCTR","21-Reference")] <- NULL
sampleNames <- c("iHT_03iCTR","iHT_90iOBS","iHT_77iOBS","aHT_1662","aHT_1838","aHT_1843","aHT_2266","iHT_02iCTR","aHT_2884","iHT_87iCTR","iHT_201iCTR","iHT_25iCTR","iHT_688iCTR","iHT_80iCTR","iHT_74iOBS","iHT_03iOBS")

### Format samples into a list
sampleTranscriptomeList <- list()
for (sampleNumber in 1:length(TPMdata)) {
  currentDataframe <- TPMdata[sampleNumber]
  order <- order(currentDataframe[1], decreasing=TRUE)
  sampleTranscriptomeList[[sampleNumber]] <- data.frame(currentDataframe[order,],row.names=row.names(currentDataframe)[order])
}
names(sampleTranscriptomeList) <- sampleNames

### Add references to a list
references <- convertIDs(references)
referenceTranscriptomeList <- list()
for (tissueNumber in 2:length(references)) {
  currentDataframe <- references[tissueNumber]
  order <- order(currentDataframe[1],decreasing = TRUE)
  referenceTranscriptomeList[[(tissueNumber-1)]] <- data.frame(currentDataframe[order,],row.names=row.names(currentDataframe)[order])
}
names(referenceTranscriptomeList) <- referenceNames

########################################################################
### Join sample list and reference list
########################################################################

transcriptomeList <- append(referenceTranscriptomeList,sampleTranscriptomeList)

########################################################################
### Subsample sample number
########################################################################

transcriptomeList <- transcriptomeList[-c(1,2,4,5,7,22,24,25,26,27,28,29,30,31,34,35,36,43,48,49,51,53)] # brain, skin, fibroblast
#transcriptomeList <- transcriptomeList[-c(22,49,53)] # no Testes blood, lympocyte
names(transcriptomeList)

########################################################################
### Subsample length
########################################################################

specifiedEnd <- 1000
for (elementNumber in 1:length(transcriptomeList)) {
  currentDatafame <- transcriptomeList[[elementNumber]]
  #selectedRows <- currentDatafame[,1] > 1 
  selectedRows <- 1:specifiedEnd
  transcriptomeList[[elementNumber]] <- data.frame(currentDatafame[1][selectedRows,],row.names=row.names(currentDatafame)[selectedRows])
}

info <- print(str(transcriptomeList))

########################################################################
### Generate list of unique ids
########################################################################

allIDs <- c()
for (df in transcriptomeList) {
  allIDs <- c(allIDs,row.names(df))
}
uniqueIDs <- unique(allIDs)
length(uniqueIDs)

### Generate a dataframe with each of the unique IDs
transcriptsDF <- data.frame(uniqueIDs)
for (df in transcriptomeList) {
  unusedIDs <- setdiff(uniqueIDs,row.names(df))
  newDF <- data.frame(c(df[,1],rep(0,length(unusedIDs))))
  row.names(newDF) <- c(row.names(df),unusedIDs)
  newDF <- newDF[order(row.names(newDF)),]
  transcriptsDF[length(transcriptsDF)+1] <- newDF
}

row.names(transcriptsDF) <- transcriptsDF[,1]  # Assign row names
transcriptsDF <- transcriptsDF[2:length(transcriptsDF)]  #Remove first column of DF
names(transcriptsDF) <- names(transcriptomeList)  # Assign column names
transcriptsMatrix <- as.matrix(transcriptsDF)  # Convert to matrix

########################################################################
### Compare, v3
########################################################################

comparisonMatrix <- rcorr(transcriptsMatrix, type="spearman")[[1]]
comparisonMatrix <- round(comparisonMatrix*100,0)

########################################################################
### Reorder, v3
########################################################################

sample.order <- comparisonMatrix["iHT_02iCTR_S13",]
sample.order.new <- names(sort(test,decreasing = TRUE))
test <- comparisonMatrix[sample.order.new,]
test <- test[,sample.order.new]

colnames(comparisonMatrix)
sample.order.num <- c(32,33,34,39,41,42,43,44,45,46,47,35,36,37,38,40,11,3,10,13,5,9,15,14,
                      25,6,7,23,1,18,19,22,2,30,16,31,20)
test <- comparisonMatrix[sample.order.num,]
comparisonMatrix <- test[,sample.order.num]

######################################################################################
### Plot heatmap
######################################################################################

start = 30
end = 100
nColors = (end - start) * 2
my_palette <- rev(colorRampPalette(c("black","red","orange","yellow","white"))(n = nColors))
my_palette <- rev(heat.colors(nColors))

col_breaks <- seq(start,end,(end-start)/(nColors))

lmat = rbind(c(0,3),c(2,1),c(0,4))
lwid = c(1,7)
lhei = c(1,5,1)

title <- paste("Spearman comparison of samples against GTEx references")

heatmap.2(comparisonMatrix,
          main = title, # heat map title
          cellnote = comparisonMatrix,
          notecol = "gray40",
          density.info="none",  # turns off density plot inside color legend
          notecex=0.7,
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette,       # use on color palette defined earlier 
          distfun=dist,
          margins =c(12,12),     # widens margins around plot
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",             # turn off column clustering
          lmat = lmat, lwid = lwid, lhei = lhei
          )
          #cellnote = comparisonMatrix,  # same data set for cell labels
          #notecol="black",      # change font color of cell labels to black



########################################################################
### Export plot
########################################################################

setwd("z:/Uthra/HT paper/Bioinformatics figures/Spearman heatmaps/")

filename <- "heatmap-10k-reordered-recolored"

png(filename=paste0(filename,".png"), 
    type="cairo",
    units="in", 
    width=14, 
    height=10, 
    pointsize=12, 
    res=400)
heatmap.2(comparisonMatrix,
          main = title, # heat map title
          cellnote = comparisonMatrix,
          notecol = "gray40",
          density.info="none",  # turns off density plot inside color legend
          notecex=0.7,
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette,       # use on color palette defined earlier 
          distfun=dist,
          margins =c(12,12),     # widens margins around plot
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          Colv="NA",             # turn off column clustering
          lmat = lmat, lwid = lwid, lhei = lhei
)
dev.off()


tiff(filename=paste0(filename,".tiff"), 
    type="cairo",
    units="in", 
    width=14, 
    height=10, 
    pointsize=12, 
    res=300)
heatmap.2(comparisonMatrix,
          main = title, # heat map title
          cellnote = comparisonMatrix,
          notecol = "gray40",
          density.info="none",  # turns off density plot inside color legend
          notecex=0.6,
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette,       # use on color palette defined earlier 
          distfun=dist,
          margins =c(12,12),     # widens margins around plot
          breaks=col_breaks,    # enable color transition at specified limits
          lmat = lmat, lwid = lwid, lhei = lhei)
dev.off()

#transcriptsDF <- addGene(transcriptsDF)

write.csv(comparisonMatrix,paste0(filename,"-comparison_matrix.csv"))
write.csv(transcriptsDF,paste0(filename,"-expression.csv"))
