### Opening external transcriptomes -- Andrew R Gross -- 2016-11-22
### Open files from various sources and organize them


########################################################################
### Header
########################################################################


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
### Import
########################################################################

### iPSC Cardiomyocites
### Burridge, Nat. Med. 2016 (Wu) | RNA-seq, GSE76314
df.cardio.bur <- read.csv("C:/Users/grossar/Bioinform/DATA/rna_seq/Reference_transcriptomes/cardio-burridge/GSE76314_Burridge_NM_Processed_Data_File_log_fpkm_plus_1_.csv", row.names = 1)
### Burridge, Nat. Med. 2016 (Wu) | Affy microarray, GSE79413

### Matsa, Cell Stem Cell, 2016 (Wu) | RNA-seq, GSE83668
df.cardio.mat <- read.delim("C:/Users/grossar/Bioinform/DATA/rna_seq/Reference_transcriptomes/cardio-matsa/GSE83668_RNASeq.txt")
df.cardio.mat <- df.cardio.mat[c(1,14,15,16,17,18)][-(1:20),]
row.names(df.cardio.mat) <- df.cardio.mat$Gene.Symbol
df.cardio.mat <- df.cardio.mat[-1]

### Kodo, Nature Cell Biology, 2016 (Wu) | RNA-seq, GSE63161
df.cardio.kod <- read.delim("C:/Users/grossar/Bioinform/DATA/rna_seq/Reference_transcriptomes/cardio-kodo/GSE63161_FPKM-all-samples.txt", stringsAsFactors = FALSE)


########################################################################
### Format
########################################################################

### Kodo, Nature Cell Biology, 2016 (Wu) | RNA-seq, GSE63161
for(line.pos in 1:nrow(df.cardio.kod)){
  line <- df.cardio.kod[line.pos,]
  ids <- strsplit(as.character(unlist(line[2])),',')[[1]]
  df.cardio.kod[line.pos,][1] <- ids[1]
  if(length(ids)>1){
    df.cardio.kod[line.pos,][2] <- ids[2]
  }
}

unique.kod <- unique(df.cardio.kod[,1])

test <- match(unique.kod,df.cardio.kod[,1])
un.df.cardio.kod <- df.cardio.kod[test,]
row.names(un.df.cardio.kod) <- un.df.cardio.kod[,1]
df.cardio.kod <- un.df.cardio.kod[-(1:2)]
########################################################################
### Sort
########################################################################

df.cardio.bur <- addMedSD(df.cardio.bur)
df.cardio.bur <- sortByMed(df.cardio.bur)
  
df.cardio.mat <- addMedSD(df.cardio.mat)
df.cardio.mat <- sortByMed(df.cardio.mat)

df.cardio.kod <- addMedSD(df.cardio.kod)
#df.cardio.kod <- data.frame(df.cardio.kod[1],addMedSD(df.cardio.kod[3:10]))
df.cardio.kod <- sortByMed(df.cardio.kod)


########################################################################
### Calculate overlaps
########################################################################

### Calculate overlaps of starting datasets
nrow(df.cardio.bur)
bur <- row.names(df.cardio.bur)
nrow(df.cardio.mat)
mat <- row.names(df.cardio.mat)
nrow(df.cardio.kod)
kod <- row.names(df.cardio.kod)

length(intersect(bur,mat))/10305
length(intersect(mat,kod))/10305
length(intersect(kod,bur))/21616

length(intersect(bur[1:5000],mat[1:5000]))/5000
length(intersect(mat[1:1000],kod[1:1000]))/1000
length(intersect(kod[1:1000],bur[1:1000]))/1000

genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters='external_gene_name', values=test.list[1:1000], mart=ensembl)

########################################################################
### Generate list of unique ids
########################################################################

TranscriptomeList <- list()
transcriptomeList <- list(df.cardio.bur,df.cardio.mat,df.cardio.kod)

for (sampleNumber in 1:length(TPMdata)) {
  currentDataframe <- TPMdata[sampleNumber]
  order <- order(currentDataframe[1], decreasing=TRUE)
  sampleTranscriptomeList[[sampleNumber]] <- data.frame(currentDataframe[order,],row.names=row.names(currentDataframe)[order])
}

########################################################################
### Generate list of unique ids
########################################################################

allIDs <- c(bur,mat,kod)
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


