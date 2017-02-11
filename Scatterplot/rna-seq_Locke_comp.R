### RNA Seq preliminary analysis, Andrew R Gross, 2016-05-16
### This script is intended to upload normalized expression data and plot a variety of features in order to make general assessments of the data

########################################################################
### Header
########################################################################

library(gplots)
library(RColorBrewer)
library(biomaRt)
library(DESeq2)
library(Hmisc)

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
########################################################################

addMedSD <- function(dataframe) {                 # A function to add a column for median value and a column for std. dev.
  median <- apply(dataframe,1,median)
  sd <- apply(dataframe,1,sd)
  return(data.frame(dataframe,median,sd))  
}
sortByMed <- function(dataframe) {                # A function to sort a dataframe by a median value
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
addDescription <- function(dataframe) {
  descriptions <- getBM(attributes=c('ensembl_gene_id','description'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  descriptions <- descriptions[match(row.names(dataframe),descriptions[,1]),]
  Descrip <- c()
  for (rowNumber in 1:length(descriptions[,1])) {
    fullDescr <- descriptions[rowNumber,][,2]
    shortDescr <- strsplit(fullDescr," \\[Source")[[1]][1]
    Descrip <- c(Descrip, shortDescr)
  }
  dataframe[length(dataframe)+1] <- Descrip
  names(dataframe)[ncol(dataframe)] <- "Descrip"
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

# Metadata
metaData <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/samples.csv", row.names=1)

# Normalized
#TPMdata <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)

# Raw, Non-normalized
#counts <- read.table("z://Data/RNAseq HT neurons and tissue/2nd rerun/20160317_Sareen_rerun-29299281.all.count", row.names=1, header=TRUE)

# Prepared DE files
importFile1 <- "z://Data/RNAseq HT neurons and tissue/2nd rerun/DE_genes_ANDREW/Obs_iHT_vs_Ctr_iHT--10000.csv"
#importFile2 <- "z://Data/RNAseq HT neurons and tissue/2nd rerun/DE_gene_list before filter/OBS-iHT_v_CTR-iHT_DE-Descr-Tues-Jul-05.csv"
#"z://Data/RNAseq HT neurons and tissue/2nd rerun/DE_gene_list after_filter/OBS_vs_CTR/OBS-iHT_vs_CTR-iHT_filtered_DE_csv_results.csv"
deOBSvCTR <- read.csv(importFile1,row.names=1,header =TRUE)

# Locke gene list
lockeGenes <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/Locke-gene-lists.csv")
locke <- lockeGenes$Ensembl_ID

########################################################################
### Formatting
########################################################################

deOBSvCTR <- deOBSvCTR[c(1,2,6,7)]                               # Remove unwanted columns
fileName <- strsplit(importFile1,"/")[[1]][length(strsplit(importFile1,"/")[[1]])]  # Declare a file name for outputting
pValueCutoff <- 0.05                               # Assign a p-value cutoff

########################################################################
### Expression assignment
########################################################################

expression = rep('re',length(deOBSvCTR[,1]))                               # Establish a vector of 're', for 'regular expression'
highP <- deOBSvCTR$p.adj < pValueCutoff                               # Generate a vector indicating which positions have p-values above the cutoff
expression[highP] <- 'de'                               # Tag all positions with differential expression 'de'.
summary(as.factor(expression))                               # Report the number of 'de' and 're' rows

deOBSvCTR[ncol(deOBSvCTR)+1] <- expression                               # Add a new column for 'expression'
names(deOBSvCTR)[ncol(deOBSvCTR)] <- "expression"                               # Name the new column

########################################################################
### Locke assignment
########################################################################

sharedGenes <- intersect(row.names(deOBSvCTR),locke)                               # Create a vector of shared genes

lockeCol <- rep("NIL",length(deOBSvCTR[,1]))                               # Create a vector of 'NIL' for 'Not in Locke'
deOBSvCTR[ncol(deOBSvCTR)+1] <- lockeCol                               # Create a new column of the 'NIL' vector
names(deOBSvCTR)[ncol(deOBSvCTR)] <- "locke"                               # Rename the new row
deOBSvCTR[sharedGenes,][ncol(deOBSvCTR)] <- "locke"                               # Change the values in the rows which were found in Locke to reflect this

deOBSvCTR$expression <- as.factor(deOBSvCTR$expression)                               # Reclassify the new columns as factors
deOBSvCTR$locke <- as.factor(deOBSvCTR$locke)                               

deOBSvCTR$OBS.median <- log10(deOBSvCTR$OBS.median+1)                               # Log transform the values
deOBSvCTR$CTR.median <- log10(deOBSvCTR$CTR.median+1)                               # 

########################################################################
### Generate specific data sets for each layer
########################################################################

deOvC_re <- deOBSvCTR[deOBSvCTR$expression=='re',]  # Declare a dataframe containing regularly expressed genes
deOvC_de <- deOBSvCTR[deOBSvCTR$expression=='de',]  # Declare a dataframe containing differentially expressed genes

deOvC_re_locke <- deOvC_re[deOvC_re$locke=='locke',]  # Declare a dataframe of regularly expressed Locke genes
deOvC_de_locke <- deOvC_de[deOvC_de$locke=='locke',]  # Declare a dataframe of differentially expressed Locke genes
deOvC_de_locke$label <- deOvC_de_locke$OBS.median

up_down_sort <- deOvC_de_locke[,2]-deOvC_de_locke[,1] # Create list of expression change to sort
deOvC_de_locke_up <- deOvC_de_locke[up_down_sort>0,]
deOvC_de_locke_down <- deOvC_de_locke[up_down_sort<0,]

deOvC_de_top10 <- deOvC_de[1:10,]   # Declare a dataframe of the top 10 most differentially expressed genes by p-value
up_down_sort <- deOvC_de[1:10,][,2]-deOvC_de[1:10,][,1]  # Sort top ten genes by direction
deOvC_de_top10_up <- deOvC_de_top10[up_down_sort>0,]  # Declare a dataframe of the top ten genes which are upregulated
deOvC_de_top10_down <- deOvC_de_top10[up_down_sort<0,]  # Declare a dataframe of the top ten genes which are down regulated.

########################################################################
### Plot
########################################################################

### Plot points
g <- ggplot() + geom_point(data=deOvC_re,aes(OBS.median,CTR.median),color="grey90",size=1.5)+theme_bw() +   # Plot the reg. expr. genes
  geom_point(data=deOvC_de,aes(OBS.median,CTR.median,size= -log10(p.adj),color=(p.adj))) +   # plot the diff. expr. genes
  geom_point(data=deOvC_re_locke,aes(OBS.median,CTR.median),color="pink",size=1.5) +   # Plot the reg. expr. Locke genes
  geom_point(data=deOvC_de_locke,aes(OBS.median,CTR.median),color="red",size=5,shape=18)   # plot the diff. expr. Locke genes

### Plot text
(h <- g + geom_text(data=deOvC_de_top10_up,aes(OBS.median,CTR.median,label=Gene),hjust=1.27,vjust=0.5,size=4) + # Plot select upreg.
  geom_text(data=deOvC_de_top10_down,aes(OBS.median,CTR.median,label=Gene),hjust=-0.28,vjust=.8,size=4) +  # Plot select downreg.
  geom_text(data=deOvC_de_locke_up,aes(x = label - 1, y = CTR.median,label=Gene),color="red",size=5) +
  geom_segment(data=deOvC_de_locke_up,aes(x = label - .2, y = CTR.median, xend = label - .8, yend = CTR.median),linetype="dotted",color="red") +
  geom_text(data=deOvC_de_locke_down,aes(x = label + 1 , y = CTR.median,label=Gene),color="red",size=5) + 
  geom_segment(data=deOvC_de_locke_down,aes(x = OBS.median + .2, y = CTR.median, xend = label + .8, yend = CTR.median),linetype="dotted",color="red"))

### Adjust placement of text
deOvC_de_locke_down
deOvC_de_locke_down$label[3] <- 3.501607 + .15

g <- h
  
### Add title and labels
g <- g + ggtitle(("Median expression of 10,000 genes\nin control vs. obese samples")) +
  scale_colour_distiller(direction = -1,"P-Value")

### Format background
g <- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

### Format margins and axis labels
h <- g + theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
               axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
               axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
               plot.margin = unit(c(1,1,1,1), "cm"),
               axis.text = element_text(size = 12)) +
              labs(x = "Expression in obese samples, log10(TPM)", y = "Expression in control samples, log10(TPM)") +
              xlim(0,5.7) + ylim(0,5.7) + scale_size_continuous(name = "P-value (-log)")

h # Preview

########################################################################
### Export plot
########################################################################

setwd("z:/Uthra/HT paper/Bioinformatics figures/Scatter plots/")

filename <- "Scatterplot_05"

png(filename=paste0(filename,".png"), 
    type="cairo",
    units="in", 
    width=14, 
    height=14, 
    pointsize=12, 
    res=300)
h
dev.off()

tiff(filename=paste0(filename,".tiff"), 
    type="cairo",
    units="in", 
    width=14, 
    height=14, 
    pointsize=12, 
    res=300)
h
dev.off()




