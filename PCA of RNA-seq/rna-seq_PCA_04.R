### RNA Seq PCA analysis, Andrew R Gross, 2016-08-30
### Input: Normalized expression data
### Output: PCA plots; tables of genes contributing to various loadings

########################################################################
### Header
########################################################################

library(gplots)
library(RColorBrewer)
library(biomaRt)
library(DESeq2)
library(Hmisc)
library(ggbiplot)
library(rgl)

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
  #print(Descrip)
  dataframe[length(dataframe)+1] <- Descrip
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
metadata.all <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/HT_plus_reference_metadata.csv",row.names=1)
#metadata.sam <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/samples.csv", row.names=1)
#metadata.ref <- read.csv("c://Users/grossar/Bioinform/DATA/rna_seq/reference_metadata.csv",row.names=1)

# Normalized
TPMdata <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)

# Raw, Non-normalized
counts <- read.table("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.all.count", row.names=1, header=TRUE)

# References
references <- read.table("c://Users/grossar/Bioinform/DATA/rna_seq/Reference_transcriptomes/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct",sep="\t",header=TRUE,row.names=1)

# Example data
#dds <- DESeqDataSet(se, design = ~ cell + dex)

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
### Formating RNA-seq data
########################################################################

### Rename rows and columns
counts <- convertIDs(counts)     # Correct row names
names(counts)
names(TPMdata)
names(counts) <- names(TPMdata)   # Correct column names
counts <- counts[row.names(metadata.sam)]

### Correct TPM data
TPMdata <- convertIDs(TPMdata)  # Remove decimals from IDs
TPMdata <- TPMdata[row.names(metadata.sam)]  # Reorder columns to match metadata

### Remove unwanted column
TPMdata <- TPMdata[-match("Reference",names(TPMdata))]
counts <- counts[-match("Reference",names(counts))]
metadata.all <- metadata.all[-match("Reference",row.names(metadata.all)),]

########################################################################
### Formating gTex reference data
########################################################################

### Remove unwanted row, correct row names
references <- references[2:length(references)]
references <- convertIDs(references)
names(references) <- referenceNames

row.names(metadata.ref) <- referenceNames
row.names(metadata.all)

### Rescale RPKM to TPM
for (columnNumber in 1:length(references)) {
  references[columnNumber] <- references[columnNumber]*1000000/sum(references[columnNumber])
}

### Round values to two decimal places
references <- round(references,2)

########################################################################
### Joining 
########################################################################

### Generate list of all ensemble IDs in each of the two sets
referenceIDs <- row.names(references)
sampleIDs <- row.names(TPMdata)

### Genereate list of the IDs not present in each
notSampleIDs <- setdiff(referenceIDs,sampleIDs)
notReferenceIDs <- setdiff(sampleIDs,referenceIDs)

### Create empty rows and append them to sample DF
sampleAddDF <- data.frame(mat.or.vec(length(notSampleIDs),length(TPMdata)))
names(sampleAddDF) <- names (TPMdata)
row.names(sampleAddDF) <- notSampleIDs
TPMdata.plus <- rbind(TPMdata,sampleAddDF)

### Create empty rows and append them to reference dataframe
referenceAddDF <- data.frame(mat.or.vec(length(notReferenceIDs),length(references)))
names(referenceAddDF) <- names(references)
row.names(referenceAddDF) <- notReferenceIDs
reference.plus <- rbind(references,referenceAddDF)

### Reorder sample DF
samplesDF <- TPMdata.plus[order(row.names(TPMdata.plus)),]
referencesDF <- reference.plus[order(row.names(reference.plus)),]

### Append references to sample DF
sam.plus.ref.df <- cbind(samplesDF,referencesDF)

########################################################################
### Filtering
########################################################################
### Metadata files are used to group PCA icons

metadata.sam <- metadata.all[!metadata.all$Source == 'Gtex',]  # Select metadata from iHT samples
metadata.sam.iHT <- metadata.all[metadata.all$Group == 'iHT',]  # Select metadata from iHT samples
metadata.sam.iHT.f <- metadata.sam.iHT[metadata.sam.iHT$Sex =='F',]  # Select female iHT samples
metadata.sam.HT <- na.omit(metadata.all[metadata.all$Type == "HT",])  # Select HT samples
metadata.sam.HT.f <- metadata.sam.HT[metadata.sam.HT$Sex =='F',]  # Select female iHT samples
metadata.sam.iPSC <- metadata.all[metadata.all$Source == "iPSC",]  # Select iPSC samples

metadata.gTex <- metadata.all[metadata.all$Source == "Gtex",]
metadata.gTex.brain <- metadata.gTex[metadata.all$Organ == "Brain",]

metadata.brain <- metadata.all[metadata.all$Organ == "Brain",]

########################################################################
### Subsetting
########################################################################

### Normalized
DF.norm.all <- sam.plus.ref.df[match(row.names(metadata.sam),names(sam.plus.ref.df))]
DF.norm.iHT <- sam.plus.ref.df[match(row.names(metadata.sam.iHT),names(sam.plus.ref.df))]
DF.norm.HT <- sam.plus.ref.df[match(row.names(metadata.sam.HT),names(sam.plus.ref.df))]
DF.norm.HT.f <- sam.plus.ref.df[match(row.names(metadata.sam.iHT.f),names(sam.plus.ref.df))]
DF.norm.iPSC <- sam.plus.ref.df[match(row.names(metadata.sam.iPSC),names(sam.plus.ref.df))]

### Raw
DF.raw <- counts[match(row.names(metadata.sam),names(counts))]
DF.raw.iHT <- counts[match(row.names(metadata.sam.iHT),names(counts))]
DF.raw.HT <- counts[match(row.names(metadata.sam.HT),names(counts))]
DF.raw.iHT.f <- counts[match(row.names(metadata.sam.iHT.f),names(counts))]
DF.raw.HT.f <- counts[match(row.names(metadata.sam.HT.f),names(counts))]
DF.raw.iPSC <- counts[]

### Refernces
DF.references <- sam.plus.ref.df[match(row.names(metadata.gTex),names(sam.plus.ref.df))]
DF.ref.brain <- sam.plus.ref.df[match(row.names(metadata.gTex.brain),names(sam.plus.ref.df))]

### Both
DF.both.brain <- sam.plus.ref.df[match(row.names(metadata.brain),names(sam.plus.ref.df))]

########################################################################
### Convert to matricies
########################################################################

### Normalized
countMatrix <- as.matrix(DF.norm.all) ; metadata.selected <- metadata.sam ; set.title <- "All samples" ; groups <- c("Type","Source","Sex","Disease") 
countMatrix <- as.matrix(DF.norm.iHT) ; metadata.selected <- metadata.sam.iHT ; set.title <- "iHT samples" ; groups <- c("Sex","Disease")
countMatrix <- as.matrix(DF.norm.HT) ; metadata.selected <- metadata.sam.HT ; set.title <- "HT samples" ; groups <- c("Source","Sex","Disease") 
countMatrix <- as.matrix(DF.norm.HT.f) ; metadata.selected <- metadata.sam.iHT.f ; set.title <- "HT samples-female" ; groups <- c("Source","Disease") 
countMatrix <- as.matrix(DF.norm.iPSC) ; metadata.selected <- metadata.sam.iPSC ; set.title <- "iPSC samples" ; groups <- c("Type","Sex","Disease") 

### Raw
countMatrix <- as.matrix(DF.raw) ; metadata.selected <- metadata.sam
countMatrix <- as.matrix(DF.raw.iHT) ; metadata.selected <- metadata.sam.iHT
countMatrix <- as.matrix(DF.raw.HT) ; metadata.selected <- metadata.sam.HT
countMatrix <- as.matrix(DF.raw.iHT.f) ; metadata.selected <- metadata.sam.iHT.f
countMatrix <- as.matrix(DF.raw.HT.f) ; metadata.selected <- metadata.sam.HT.f
countMatrix <- as.matrix(DF.raw.iPSC) ; metadata.selected <- metadata.sam.iPSC

### References
countMatrix <- as.matrix(DF.references) ; metadata.selected <- metadata.gTex ; set.title <- "All samples"
countMatrix <- as.matrix(DF.ref.brain) ; metadata.selected <- metadata.brain.gTex

### Both
countMatrix <- as.matrix(DF.both.brain) ; metadata.selected <- metadata.brain

###
countMatrix <- round(countMatrix)
head(countMatrix)
(metadata.selected <- droplevels(metadata.selected))

####################################################################################################################
######### Plotting PCA
####################################################################################################################
### DESeq data -> plotPCA --- Good, but doesn't have labels
### DESeq data -> ggplot --- Appears to be best so far.  Makes previous method obsolete
### prcomp method -> biplot --- Useless.  Ugly.
### prcomp method -> ggbiplot --- Doesn't seem to work
### prcomp method -> ggplot --- Good.  Needs formatting
### prcomp method -> plot3d --- 
### princomp -> any of these  --- Doesn't appear to be compatible with my dataset size

########################################################################
### Calculate components using DESeq
########################################################################

### Generate DESeq object from HT counts
deHT <- DESeqDataSetFromMatrix(countData = countMatrix, colData = metadata.selected, design = ~ Sex + Disease)
deHT <- DESeqDataSetFromMatrix(countData = countMatrix, colData = metadata.selected, design = ~ Source + Sex)
deHT <- DESeqDataSetFromMatrix(countData = countMatrix, colData = metadata.selected, design = ~ Organ)

nrow(deHT)                                  # Count rows in deHT
deHT <- deHT[rowSums(counts(deHT)) > 1, ]   # Retain only rows with at least one count
nrow(deHT)                                  # Count rows again
#deHTtrans <- rlog(deHT,blind=FALSE)         # Perform rlog transformation
deHTtrans <- log(deHT)  # Perform simple log transformation

setwd("c:\\Users/grossar/Bioinform/DATA/rna_seq/")
filename <- "PCA_deseq_gtex.txt"
#write.csv(deHTtrans,file=filename)
#writebin

### Collect plot details
(deseq.pca.data <- plotPCA(deHTtrans, intgroup = c("Source", "Disease", "Sex"), returnData=TRUE))
(deseq.pca.data <- plotPCA(deHTtrans, intgroup = "Organ", returnData=TRUE))

percentVar <- round(100 * attr(deseq.pca.data, "percentVar"))

######################################################################################
### Calculate components using prcomp
######################################################################################

### Filter out empty rows
nrow(countMatrix)                                  # Count rows
countMatrix <- countMatrix[rowSums(countMatrix) > 1, ]   # Retain only rows with at least one count
nrow(countMatrix)                                  # Count rows again

### Perform rlog transformation and transpose
countMatrix.log <- rlog(countMatrix)  # Perform log transformation
count.trans <- t(countMatrix.log)                 # Transpose matrix to place samples as rows and genes as columns

### Perform PCA with prcomp
count.pca <- prcomp(count.trans, center = TRUE, scale. = TRUE) # scale. = TRUE is highly advisable, but default is FALSE.
#count.pca <- prcomp(t(assay(deHTtrans)), center = TRUE, scale. = TRUE)   # prcomp using the DESeq data

### Display component info
plot(count.pca, type = "l")  # A plot of the size of the components
summary(count.pca)  # An out put of summarized statistics
#print(count.pca)   # Components.  Too large to view
#predict(ir.pca,  newdata=tail(log.ir, 2))  # I don't know what this does. It's written for iris data anyway.

######################################################################################
### Visualize transformation
######################################################################################

### Examine transformation effect
#par( mfrow = c( 1, 2 ) )

#plot(assay(deHTtrans)[,1:2], pch=16, cex=0.3)
#plot(log2(countMatrix[,1:2]), pch=16, cex=0.3)     # Plot the log transformed values
#plot(countMatrix.log[,1:2], pch=16, cex=0.3)           # Plot the log transformed values


######################################################################################
### Plot DESeq data using ggplot
######################################################################################

### Plot deseq.pca.data using ggplot
ggplot(deseq.pca.data, aes(PC1, PC2, color=Source, shape=Disease)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=row.names(deseq.pca.data)),hjust=0,vjust=0)

ggplot(deseq.pca.data, aes(PC1, PC2, color=Organ)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=row.names(deseq.pca.data)),hjust=0,vjust=0)

######################################################################################
### Plot prcomp data using ggplot
######################################################################################

### Specify components
count.pca.df <- data.frame(count.pca$x)
count.pca.df <- deseq.pca.data

ggplot(count.pca.df, aes(PC1, PC2)) + geom_point(size=3) +
  geom_text(aes(label=row.names(count.pca.df)),hjust=0,vjust=0)

ggplot(count.pca.df, aes(PC1, PC2, color=deseq.pca.data$Source, shape=deseq.pca.data$Disease)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=row.names(deseq.pca.data)),hjust=0,vjust=0)

ggplot(count.pca.df, aes(PC1, PC2, color=deseq.pca.data$Organ)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=row.names(deseq.pca.data)),hjust=0,vjust=0)

######################################################################################
### Plot prcomp data in 3D
######################################################################################

plot3d(count.pca.df[,1:3])
text3d(count.pca.df[,1:3],texts=rownames(count.trans))






######################################################################################
### Calculate loading vectors
######################################################################################









######################################################################################
### Examples, scratchwork
######################################################################################

### Calculate vectors for PC
coords <- NULL
for (i in 1:nrow(pc$loadings)) {
  coords <- rbind(coords, rbind(c(0,0,0),pc$loadings[i,1:3]))
}

coords <- NULL
for (i in 1:nrow(pc2$rotation)) {
  coords <- rbind(coords, rbind(c(0,0,0),pc$rotation[i,1:3]))
}

# Plot points, labels, component labels, and component vectors
plot3d(pc$scores[,1:3])
text3d(pc$scores[,1:3],texts=rownames(test.trans))
text3d(pc$loadings[,1:3], texts=rownames(pc$loadings), col="red")
lines3d(coords, col="red", lwd=4)


### some kind of validation step
#plotMA(res, main="DESeq2", ylim=c(-2,2))
#plotCounts(deHT, gene=which.min(res$padj), intgroup="Disease")



########################################################################
### rLog transforming dataset
########################################################################

sam.plus.ref.df.rlog <- rlog(round(as.matrix(sam.plus.ref.df)))

sam.plus.ref.df.log <- log(as.matrix(sam.plus.ref.df))


#setwd("c:\\Users/grossar/Bioinform/DATA/rna_seq/")
filename <- "samples-plus-ref_rlog-transformed"
#write.csv(sam.plus.ref.df.rlog,file = paste0(filename,".csv"))
#sam.plus.ref.df.rlog <- read.csv(sam.plus.ref.df.rlog,file = paste0(filename,date,".csv"))


########################################################################
### TEMPORARY -- vALIDATE
########################################################################

### subset, then transform.  Compare to transform, then subset.

test <- round(as.matrix(sam.plus.ref.df[c(1,10,20,30)]))
test.rlog <- rlog(test)

test.compare <- sam.plus.ref.df.rlog[,c(1,10,20,30)]

### Examine transformation effect
par( mfrow = c( 1, 2 ) )
plot(test.rlog[,1:2], pch=16, cex=0.3)
plot(test.compare[,1:2], pch=16, cex=0.3)

test.log <- log(test)
test.compare <- sam.plus.ref.df.log[,c(1,10,20,30)]

plot(test.log[,1:2], pch=16, cex=0.3)
plot(test.compare[,1:2], pch=16, cex=0.3)

