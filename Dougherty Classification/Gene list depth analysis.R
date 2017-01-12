### Line graphs for depth converage -- Andrew R Gross -- 2016-11-15
### Compare the number of genes detected across a range of gene list lengths (search depth)


########################################################################
### Header
########################################################################

ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

library(ggplot2)

########################################################################
### Functions
########################################################################
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
### Data input
########################################################################

### pSI table for gtex-full 
psi.gtex.full <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/specificity.indexes/pSI_gTEX_full_10-31-16.csv", header=T,row.names=1)
HT <- psi.gtex.full[16]
title <- "GTEx, Full"
counts <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/specificity.indexes/pSI_gtex_counts_10_03_16.csv")

### pSI table for gtex, brains
psi.gtex.brain <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/specificity.indexes/psi_gtex_brain_SPEC-INDEX_11-15-16.csv", header = T, row.names = 1)
HT <- psi.gtex.brain[10]
title <- "GTEx, Brain"
counts <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/specificity.indexes/psi_gtex_brain_COUNTS_11-15-16.csv", header = T, row.names = 1)

gHT <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/gene lists/GTEx.HT-2000_MonOct311657.txt", header = FALSE)[,1]

aHT <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/gene lists/aHT-2000_MonOct311656.txt", header = FALSE)[,1]

aiHT<- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/gene lists/HT-2000_MonOct311658.txt", header = FALSE)[,1]

iHT <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/gene lists/iHT-2000_MonOct311659.txt", header = FALSE)[,1]

iMN<- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/gene lists/MN-2000_MonOct311704.txt", header = FALSE)[,1]

########################################################################
### Format
########################################################################

names(HT)[1] <- "Hypothalamus"
HT <- subset(HT, Hypothalamus <0.05)
HT <- addGene(HT)

psi.gtex.ht.05 <- subset(HT, Hypothalamus <0.05)[,2]
psi.gtex.ht.01 <-  subset(HT, Hypothalamus <0.01)[,2]
psi.gtex.ht.001 <-  subset(HT, Hypothalamus <0.001)[,2]
psi.gtex.ht.0001 <-  subset(HT, Hypothalamus <0.0001)[,2]

########################################################################
### Calculate number of gene matches
########################################################################

psi.gtex.ht <- psi.gtex.ht.0001
depths <- seq(1:500)
aHT.depths <- iHT.depths <- aiHT.depths <- gHT.depths <- iMN.depths <- c()

for (depth in depths) {
  aHT.depths <- c(aHT.depths,length(intersect(aHT[1:depth],psi.gtex.ht))/depth)
}
for (depth in depths) {
  iHT.depths <- c(iHT.depths,length(intersect(iHT[1:depth],psi.gtex.ht))/depth)
}
for (depth in depths) {
  aiHT.depths <- c(aiHT.depths,length(intersect(aiHT[1:depth],psi.gtex.ht))/depth)
}
for (depth in depths) {
  gHT.depths <- c(gHT.depths,length(intersect(gHT[1:depth],psi.gtex.ht))/depth)
}
for (depth in depths) {
  iMN.depths <- c(iMN.depths,length(intersect(iMN[1:depth],psi.gtex.ht))/depth)
}

depth.df <- data.frame(depths,aHT.depths,iHT.depths,aiHT.depths,gHT.depths,iMN.depths)


#getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters='ensembl_gene_id', values=row.names(subset(HT, Hypothalamus <0.001)), mart=ensembl)
#gene.abrv.df <- getBM(attributes=c('ensembl_gene_id','external_gene_name','hgnc_symbol'), filters='ensembl_gene_id', values=row.names(subset(HT, Hypothalamus <0.001)), mart=ensembl)

########################################################################
### Plot
########################################################################

### Plot 
ggplot(data = depth.df, aes(x = depths)) +
  geom_line(aes(y = gHT.depths), color = 'black', size = 1) +
  geom_line(aes(y = aHT.depths), color = 'blue', size = 1) + 
  geom_line(aes(y = aiHT.depths), color = 'cyan3', size = 1) +
  geom_line(aes(y = iHT.depths), color = 'darkgreen', size = 1) +
  geom_line(aes(y = iMN.depths), color = 'red', size = 1) +
  theme_light() +
  theme(plot.title = element_text(size = rel(3)), axis.text = element_text(size = rel(2))) +
  ylim(0,0.3) +
  labs(title = paste(title, ", p = 0.0001"))


