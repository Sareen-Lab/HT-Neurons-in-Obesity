### RNA Seq Tissue comparison by Spearman Correlation, Andrew R Gross, 2016-05-16
### Input: Normalized RNA-seq data; Reference expression data; metadata
### Output: Heatmap displaying pairwise similarity as calculated by Spearman correlation
###

########################################################################
### Header
########################################################################

library(gplots)
library(RColorBrewer)
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

addMedSD <- function(dataframe) {
  median <- apply(dataframe,1,median)
  sd <- apply(dataframe,1,sd)
  return(data.frame(dataframe,median,sd))  
}

sortByMed <- function(dataframe) {
  order <- order(dataframe$median,decreasing=TRUE)
  return(dataframe[order,])
}

convertIDs <- function(dataframe) {
  ensemblIDs <- c()
  for (rowName in row.names(dataframe)) {
    ensemblID <- strsplit(rowName,"\\.")[[1]][1]
    ensemblIDs <- c(ensemblIDs,ensemblID)
  }
  row.names(dataframe)<-ensemblIDs
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
### Import Data
########################################################################

### Normalized
TPMdata <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)
metaData <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/samples.csv", row.names=1)

### Import references

references <- read.table("c://Users/grossar/Bioinform/DATA/rna_seq/Reference_transcriptomes/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct",sep="\t",header=TRUE,row.names=1)

referenceNames <- c("Adipose--subcutaneous","Adipose--omentum","Adrenal_gland","Aorta","Coronary_artery",
                    "Tibial_artery","Bladder","Brain--amygdala","Brain--anteriror_cingulate","Brain--caudate_nucleous",
                    "Brain--cerebellar_hemisphere","Brain--cerebellum","Brain--cortex","Brain--Frontal_cortex",
                    "Brain--hippocampus","Brain--hypothalamus","Brain--nucleus_accumbens","Brain--putamen",
                    "Brain--spinal_cord","Brain--Substantia_nigra","Mammary","Lymphocyte","Fibroblast","Ectocervix",
                    "Endocervix","Colon--sigmoid","Colon--transverse","Gastroesophageal_junction","Esophagus--mucosa",
                    "Esophagus--muscularis","Fallopian_tube","Heart--Atrial","Heart--left_ventricle",
                    "Kidney","Liver","Lung","Salvitory_gland","Skeletal_muscle","Tibial_nerve","Ovary","Pancreas",
                    "Pituitary","Prostate","Skin--Suprapubic","Skin--Leg","Small_intestine","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_blood")

########################################################################
### Formating
########################################################################

### Replace row and column names
sampleNames <- c("iHT_03iCTR","iHT_90iOBS","iHT_77iOBS","iHT_02iOBS","iMN_87iCTR","iMN_201iCTR","aHT_1662_S6","aHT_1838_S7","aHT_1843_S8","aHT_2266_S9","iHT_02iCTR_S13","aHT_2884_S11","21-Reference","iHT_87iCTR","iHT_201iCTR","iHT_25iCTR_S16","iHT_688iCTR_","iHT_80iCTR","iHT_74iOBS","iHT_03iOBS")
names(TPMdata) <- sampleNames
TPMdata <- convertIDs(TPMdata)

### Remove 02iOBS and both iMN
TPMdata[c("iHT_02iOBS","iMN_87iCTR","iMN_201iCTR","21-Reference")] <- NULL
sampleNames <- names(TPMdata)

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
### Add samples to full data list
########################################################################

transcriptomeList <- append(referenceTranscriptomeList,sampleTranscriptomeList)
#transcriptomeList <- referenceTranscriptomeList
#transcriptomeList <- sampleTranscriptomeList

########################################################################
### Subsample sample number
########################################################################

length(transcriptomeList)
transcriptomeList <- transcriptomeList[-c(22,49,53)] # no Testes blood, lympocyte
transcriptomeList <- transcriptomeList[-c(1,2,4,5,7,22,24,25,26,27,28,29,30,31,34,35,36,43,48,49,51,53)] # brain, skin, fibroblast
transcriptomeList <- transcriptomeList[[]] # Assorted
str(transcriptomeList)

########################################################################
### Subsample length
########################################################################

specifiedEnd <- 10000
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

row.names(transcriptsDF) <- transcriptsDF[,1]
transcriptsDF <- transcriptsDF[2:length(transcriptsDF)]
names(transcriptsDF) <- names(transcriptomeList)

transcriptsMatrix <- as.matrix(transcriptsDF)

########################################################################
### Compare, v3
########################################################################

comparisonMatrix <- rcorr(transcriptsMatrix, type="spearman")[[1]]
comparisonMatrix <- round(comparisonMatrix*100,0)

######################################################################################
### Plot heatmap
######################################################################################

nColors=99
start = 40
end = 100
my_palette <- colorRampPalette(c("yellow","white","blue"))(n = nColors)
my_palette <- colorRampPalette(c("red","black","green"))(n = nColors)
my_palette <- colorRampPalette(c("black","red","orange","yellow","white"))(n = nColors)
my_palette <- colorRampPalette(c("black","red","yellow","white"))(n = nColors)

col_breaks <- seq(start,end,(end-start)/(nColors))

lmat = rbind(c(0,3),c(2,1),c(0,4))
lwid = c(1,7)
lhei = c(1,5,1)

title <- paste("Spearman comparison",specifiedEnd)

heatmap.2(comparisonMatrix,
          main = title, # heat map title
          cellnote = comparisonMatrix,
          notecol = "gray40",
          density.info="none",  # turns off density plot inside color legend
          notecex=0.6,
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette,       # use on color palette defined earlier 
          distfun=dist,
          #dendrogram="column",     # only draw a row dendrogram
          #Rowv="NA",
          margins =c(12,12),     # widens margins around plot
          breaks=col_breaks,    # enable color transition at specified limits
          
          #cellnote = comparisonMatrix,  # same data set for cell labels
          #notecol="black",      # change font color of cell labels to black
          #Colv="NA"             # turn off column clustering
          lmat = lmat, lwid = lwid, lhei = lhei
          )



########################################################################
### Export plot
########################################################################

setwd("z:/Uthra/HT paper/Bioinformatics figures/Spearman heatmaps/")

filename <- "heatmap-10k-select"

png(filename=paste0(filename,".png"), 
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
