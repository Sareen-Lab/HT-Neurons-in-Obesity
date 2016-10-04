### Tissues Specific Expression Analysis, Andrew R Gross, 2016-08-26
### Using a package developed by the Dougherty lab, calculate the uniqueness of a set of genes and the probability
### of another set containing those genes.

########################################################################
### Header
########################################################################
library(biomaRt)
library(reshape2)
library(gdata)
library(ggplot2)
library(pSI)
library(RColorBrewer)
library(Hmisc)
library(grid)
library(gridExtra)
library(scales)
load("C:/Users/grossar/Bioinform/DATA/rna_seq/Dougherty/pSI.data/data/human.rda")
help(pSI)

ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

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
reorder_size <- function(x) {
  factor(x, levels = names(sort(table(x))))
}

generate.gtex.data <- function(ids.of.interest.full,length,output.table) {
  ### Generate list of ids ##########################################
  ids.of.interest <- ids.of.interest.full[1:length] # Subselects based on the specified depth
  ### Calculate p-values ############################################
  gtex.results <- fisher.iteration(pSIs = psi.gtex, candidate.genes = ids.of.interest)
  ### Output data table #############################################
  print(output.table) # Reminds user whether they requested an output table of p-values
  if(output.table == TRUE) {
    print('Outputting table')
    setwd("z://Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/pSI_results/")
    write.csv(gtex.results, paste0("gtex_",geneset,"_",length,"_",strftime(Sys.time(),"%a%b%d%H%M"),".csv"))
  }
  ### Format gtex data ##############################################
  gtex.plotting <- gtex.results[gtex.results$`0.01 - adjusted` != 1,]
  names(gtex.plotting) <- c('p0.05','p0.01','p0.001','p1e-4')
  gtex.plotting <- gtex.plotting[c('p0.01','p0.001','p1e-4')]
  gtex.plotting <- -log10(gtex.plotting)
  gtex.plotting$tissue <- row.names(gtex.plotting)
  gtex.plotting[order(gtex.plotting$p0.01,decreasing = TRUE),]
  #gtex.plotting.melt <- melt(gtex.plotting,id.var = "tissue")
  return(gtex.plotting)
}

generate.gtex.fig.from.data <- function(gtex.plotting.melt) {
  ### Plot gtex data ################################################
  gtex.plot <- ggplot(data = gtex.plotting.melt,aes(x = reorder(tissue, -value), y = value, fill = variable)) +
    geom_bar(position = "dodge",stat="identity") +
    labs(title = paste("Overlap of the top",length,"genes in",geneset,"and gtex"),
         x = "Tissue",
         y = "Negative log10 p-value") +
    theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
          axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
          axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
          plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_brewer("p-value",palette = 'Reds')
  ### Display plot ##################################################
  return(gtex.plot)
  ### Save plot #####################################################
  tiff(filename=paste0("gtex_",geneset,"_",length,"_",strftime(Sys.time(),"%a%b%d%H%M"),".tiff"), 
       type="cairo",
       units="in", 
       width=14, 
       height=14, 
       pointsize=12,
       res=100)
  gtex.plot
  dev.off() 
}

generate.tsea.figure <- function(genes.of.interest.full,length,output.table) {
  ### Generate list of ids ##########################################
  genes.of.interest <- genes.of.interest.full[1:length]
  ### Calculate p-values ############################################
  tsea.results <- fisher.iteration(pSIs = psi.tsea, candidate.genes = genes.of.interest)
  ### Output data table #############################################
  print(output.table)
  if(output.table == TRUE) {
    print('Outputting table')
    setwd("z://Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/pSI_results/")
    write.csv(tsea.results, paste0("tsea_",geneset,"_",length,"_",strftime(Sys.time(),"%a%b%d%H%M"),".csv"))
  }
  ### Format tsea data ##############################################
  tsea.plotting <- tsea.results[tsea.results$`0.01 - adjusted` != 1,]
  names(tsea.plotting) <- c('p0.05','p0.01','p0.001','p1e-4')
  tsea.plotting <- tsea.plotting[c('p0.01','p0.001','p1e-4')]
  tsea.plotting <- -log10(tsea.plotting)
  tsea.plotting$tissue <- row.names(tsea.plotting)
  tsea.plotting[order(tsea.plotting$p0.01,decreasing = TRUE),]
  tsea.plotting.melt <- melt(tsea.plotting,id.var = "tissue")
  ### Plot tsea data ################################################
  tsea.plot <- ggplot(data = tsea.plotting.melt,aes(x = reorder(tissue, -value), y = value, fill = variable)) +
    geom_bar(position = "dodge",stat="identity") +
    labs(title = paste("Overlap of the top",length,"genes in",geneset,"and tsea"),
         x = "Tissue",
         y = "Negative log10 p-value") +
    theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
          axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
          axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
          plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_brewer("p-value",palette = 'Reds')
  ### Display plot ##################################################
  return(tsea.plot)
}

add.bands.between <- function(data.frame,width) {
  new.df <- data.frame[1,]
  new.df <- rbind(new.df,rep(width,length(data.frame)))
  for(row in 2:length(data.frame[,1])) {
    new.df <- rbind(new.df,data.frame[row,])
    new.df <- rbind(new.df,rep(width,length(data.frame)))
  }
  return(new.df)
}
########################################################################
### Data input
########################################################################

# Normalized refseq data in units of TPM
TPMdata <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)

# Metadata
metadata.sam <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/samples.csv", row.names=1)

### Gtex references
references <- read.table("c://Users/grossar/Bioinform/DATA/rna_seq/Reference_transcriptomes/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct",sep="\t",header=TRUE,row.names=1)

########################################################################
### Formatting
########################################################################
### Format TPM
TPMdata <- TPMdata[row.names(metadata.sam)]  # Reorder columns to match metadata

### Remove unwanted column
TPMdata <- TPMdata[-match("Reference",names(TPMdata))]
metadata.sam <- metadata.sam[-match("Reference",row.names(metadata.sam)),]

########################################################################
### pSI determination
########################################################################

#pSI.input <- sample.data$pSI.input
#pSI.input <- references[2:ncol(references)]
#pSI.output <- specificity.index(pSI.input)
#pSI.thresholds <- pSI.list(pSI.output)
psi.counts <- pSI.count(psi.gtex)

########################################################################
### Save pSI tables
########################################################################

#setwd("z://Data/RNAseq HT neurons and tissue/Andrews_files/")
#write.csv(pSI.output,"pSI_gTEX_full_8-26-16.csv")
#write.csv(psi.counts, "z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_gTEX_counts_9-29-16.csv")

########################################################################
### Load pSI tables
########################################################################

### Dougherty tissue list, NAR 2015, Table S3 -- 25 tissues, 18056 genes labeled with gene names
psi.tsea<-read.table("http://genetics.wustl.edu/jdlab/files/2015/10/TableS3_NAR_Dougherty_Tissue_gene_pSI_v3.txt",header=T,row.names=1)

### Dougherty cell types, JoN 2016 -- 60 cell types, 16,866 genes in gene name format
psi.csea <- human$developmental.periods$psi.out

### Gtex tissues -- 53 tissues, 56,418 genes labeled with ENSB IDs
psi.gtex <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_gTEX_full_8-26-16.csv", header=T,row.names=1)
psi.gtex <- convertIDs(psi.gtex)
psi.counts <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_gTEX_counts_9-29-16.csv",row.names = 1)

########################################################################
### Subsetting
########################################################################

metadata.sam.iHT <- metadata.sam[metadata.sam$Group == 'iHT',]  # Select metadata from iHT samples
metadata.sam.iHT.f <- metadata.sam.iHT[metadata.sam.iHT$Sex =='F',]  # Select female iHT samples
metadata.sam.HT <- na.omit(metadata.sam[metadata.sam$Type == "HT",])  # Select HT samples
metadata.sam.HT.f <- metadata.sam.HT[metadata.sam.HT$Sex =='F',]  # Select female iHT samples
metadata.sam.iPSC <- metadata.sam[metadata.sam$Source == "iPSC",]  # Select iPSC samples
metadata.sam.aHT <- metadata.sam[metadata.sam$Source == "Adult",]  # Select adult HT samples
metadata.sam.MN <- na.omit(metadata.sam[metadata.sam$Type == "MN",])  # Select HT samples

HT.data <- TPMdata[match(row.names(metadata.sam.HT),names(TPMdata))]
iHT.data <- TPMdata[match(row.names(metadata.sam.iHT),names(TPMdata))]
aHT.data <- TPMdata[match(row.names(metadata.sam.aHT),names(TPMdata))]
iPSC.data <- TPMdata[match(row.names(metadata.sam.iPSC),names(TPMdata))]
MN.data <- TPMdata[match(row.names(metadata.sam.MN),names(TPMdata))]
iHT.f.data <- TPMdata[match(row.names(metadata.sam.iHT.f),names(TPMdata))]
HT.f.data <- TPMdata[match(row.names(metadata.sam.HT.f),names(TPMdata))]

datalist <- list(HT.data,iHT.data,aHT.data,iPSC.data,MN.data,iHT.f.data,HT.f.data)
names(datalist) <- c("HT","iHT","aHT","iPSC","MN","iHT.f","HT.f")

########################################################################
### Reorder genes by expression
########################################################################
### Add Median and sort by median
for(list.element in 1:length(datalist)) {
  datalist[[list.element]] <- sortByMed(addMedSD(datalist[[list.element]]))[1:4000,]
}

### Generate full lists of genes
genelist <- datalist
for(list.element in 1:length(datalist)) {
  current.genelist <- addGene(genelist[[list.element]])
  genelist[[list.element]] <- current.genelist[ncol(current.genelist)]
}

########################################################################
### Subset gene list
########################################################################

names(genelist)  # Lists the available data sets in the list of datasets
geneset <- 3  # Specifies the data set to use
genes.of.interest.full <- genelist[[geneset]][,1] # Declares the genes to use
ids.of.interest.full <- row.names(genelist[[geneset]]) # Declares the IDs of the genes to use
geneset <- names(genelist[geneset]) # Names the geneset for use in plot titles
print(geneset) # Reports the geneset in use

########################################################################
### Generate gtex plot data
########################################################################

length <- 800 

ids.of.interest <- ids.of.interest.full[1:length] # Subselects based on the specified depth
### Calculate p-values ############################################
gtex.results <- fisher.iteration(pSIs = psi.gtex, candidate.genes = ids.of.interest)
#setwd("z://Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/pSI_results/")
#write.csv(gtex.results, paste0("gtex_",geneset,"_",length,"_",strftime(Sys.time(),"%a%b%d%H%M"),".csv"))

### Format gtex data ##############################################
absent.samples <- row.names(gtex.results[gtex.results$`0.05 - adjusted` == 1,])
names(gtex.results) <- c('p0.05','p0.01','p0.001','p1e-4')
#gtex.plotting <- gtex.plotting[c('p0.01','p0.001','p1e-4')]
#gtex.plotting <- -log10(gtex.plotting)
#gtex.plotting$tissue <- row.names(gtex.plotting)
#gtex.plotting[order(gtex.plotting$p0.01,decreasing = TRUE),]

########################################################################
### Format counts and results
########################################################################

### Remove empty rows
gtex.counts.trim <- psi.counts[gtex.results$p0.05 != 1]
gtex.results.trim <- data.frame(t(gtex.results[gtex.results$p0.05 != 1,]))

### Reverse row order
#gtex.counts.trim <- gtex.counts.trim[c(4,3,2,1),]
gtex.results.trim <- gtex.results.trim[c(4,3,2,1),]

### Add edge bands
gtex.counts.trim <- add.bands.between(gtex.counts.trim,1)
gtex.counts.trim <- add.bands.between(gtex.counts.trim,0)
#gtex.results.trim <- add.bands.between(gtex.results.trim,31)
gtex.results.trim <- add.bands.between(gtex.results.trim,10^-60)

### Melt and join
gtex.plotting <- melt(gtex.counts.trim)
gtex.plotting[3] <- melt(gtex.results.trim)[2]
gtex.plotting[2] <- log10(gtex.plotting[2]+1)
gtex.plotting[3] <- -log10(gtex.plotting[3])
names(gtex.plotting) <- c('Sample', 'Set.size', 'p.Value')

########################################################################
### Plot one bullseye
########################################################################

sample.pos <- 14

seq(1,length(gtex.plotting[,1]),8)[sample.pos]
gtex.plotting.sub <- gtex.plotting[seq(1,length(gtex.plotting[,1]),8)[sample.pos]:seq(1,length(gtex.plotting[,1]),8)[sample.pos]+8,]

gtex.plotting.sub <- gtex.plotting[(((sample.pos-1)*8)+1):(((sample.pos-1)*8)+8),]
#myPalette <- colorRampPalette(c("black","white","yellow","yellow",'orange','red','firebrick4'))
myPalette <- colorRampPalette(c("gray","yellow",'orange','red','firebrick4','black'))
#myPalette <- colorRampPalette(c("black","white","yellow","yellow",'orange','orange','red','red','red4','red4'))

ggplot(data = gtex.plotting.sub, aes(x = Sample, y = Set.size, fill = p.Value)) +
  geom_bar(width = 1, stat = "identity") + 
  scale_y_continuous(limits = c(0,12)) +
  xlab(gtex.plotting.sub[1][1,]) +
  theme(axis.title.x = element_text(face="bold", size=14,margin =margin(0,0,10,0)),
        axis.title.y = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"), axis.text = element_text(size = 12),
        panel.background = element_rect(fill = 'white')) +  coord_polar() +
  scale_fill_gradientn(colors = myPalette(100), limits=c(0, 50), oob = squish)
  

########################################################################
### Experimenting with layers
########################################################################

seq(1,length(gtex.plotting[,1]),4)[sample.pos]
gtex.plotting.sub <- gtex.plotting[seq(1,length(gtex.plotting[,1]),4)[sample.pos]:seq(1,length(gtex.plotting[,1]),4)[sample.pos]+4,]

gtex.plotting.sub <- gtex.plotting[(((sample.pos-1)*4)+1):(((sample.pos-1)*4)+4),]


base <- ggplot(data = gtex.plotting.sub, aes(x = Sample, y = Set.size, fill = p.Value)) +
  geom_bar(width = 1, stat = "identity") + scale_fill_gradientn(colors = myPalette(100), limits=c(0, 50), oob = squish)

layer1 <- ggplot(data = plotting.lines, aes(x = Sample, y = Set.size)) +
  geom_bar(aes(y = Set.size/.9),width = 0.5, stat = "identity", fill = c(alpha("black",0.05),'black',alpha("black",0.5),'black',alpha("black",0.15),'black',alpha("black",0),'black'))

base <- ggplot(data = gtex.plotting.sub, aes(x = Sample, y = Set.size, fill = p.Value)) +
  geom_bar(width = 1, stat = "identity") + scale_fill_gradientn(colors = myPalette(100), limits=c(0, 50), oob = squish)

layer1 <- ggplot(data = plotting.lines, aes(x = Sample, y = Set.size)) +
  geom_bar(aes(y = Set.size/.9),width = 0.5, stat = "identity", fill = c(alpha("black",0.05),'black',alpha("black",0.5),'black',alpha("black",0.15),'black',alpha("black",0),'black'))


colour=alpha("black",0.15)

########################################################################
### Plot faceted plots
########################################################################

generate.bullseyes <- function(dataframe) {
  plot.list <- list()
  myPalette <- colorRampPalette(c("white","yellow",'orange','red','firebrick4','black'))
  #myPalette <- colorRampPalette(c("black","white","yellow","yellow",'orange','orange','red','red','red4','red4'))
  for(sample.pos in seq(1,length(dataframe[,1])/8)){
    gtex.plotting.sub <- gtex.plotting[(((sample.pos-1)*8)+1):(((sample.pos-1)*8)+8),]
    #gtex.plotting.sub <- gtex.plotting[start.row:(start.row+8),]
    current.plot <- ggplot(data = gtex.plotting.sub, aes(x = Sample, y = Set.size, fill = p.Value)) +
      geom_bar(width = 1, stat = "identity") + 
      scale_y_continuous(limits = c(0,18)) +
      xlab(gtex.plotting.sub[1][1,]) +
      theme(axis.title.x = element_text(face="bold", size=14,margin =margin(0,0,15,0)),
            axis.title.y = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
            plot.margin = unit(c(0,0,0,0), "cm"), axis.text = element_text(size = 12),
            panel.background = element_rect(fill = 'white'), legend.position="none") +  
      coord_polar() +
      scale_fill_gradientn(colors = myPalette(100), limits=c(0, 50), oob = squish)
    plot.list[[length(plot.list)+1]] <- current.plot
  }
  print(paste("List contains",length(plot.list),"plots"))
  return(plot.list)
}

plot.list <- generate.bullseyes(gtex.plotting)

########################################################################
### Select data to plot
########################################################################

p.l <- plot.list[1:15]
p.l <- plot.list[14:28]
  
(tiled.figure<-grid.arrange(p.l[[1]],p.l[[2]],p.l[[3]],p.l[[4]],p.l[[5]],
                            p.l[[6]],p.l[[7]],p.l[[8]],p.l[[9]],p.l[[10]],
                            p.l[[11]],p.l[[12]],p.l[[13]],p.l[[14]],p.l[[15]],
                            ncol = 5, top = title))


(tiled.figure <- grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],
                              ncol = 4, top = title))
nColors=99
start = 40
end = 100
(n = nColors)


  scale_colour_manual(breaks = c("0", "1", "3", "6", "9", "12"),
                      labels = c("0 month", "1 month", "3 months",
                                 "6 months", "9 months", "12 months"),
                      values = c("#E69F00", "#56B4E9", "#009E73", 
                                 "#F0E442", "#0072B2", "#D55E00"))

  scale_fill_manual(values=c("red", "blue", "green"))
 +
  


gtex.plot <- ggplot(data = gtex.plotting.melt,aes(x = reorder(tissue, -value), y = value, fill = variable)) +
  geom_bar(position = "dodge",stat="identity") + 
  labs(title = paste("Overlap of the top",length,"genes in",geneset,"and gtex"),
       x = "Tissue",
       y = "Negative log10 p-value") +
  theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer("p-value",palette = 'Reds')

gtex.plot











counts.with.bands <- add.bands.between(psi.counts,9)

#gtex.counts <- melt(psi.counts)
gtex.counts <- melt(counts.with.bands)
gtex.counts[2] <- log10(gtex.counts[2]+1)

gtex.with.bands <- add.bands.between(data.frame(t(gtex.results)),20)
temp <- melt(t(gtex.results))

gtex.counts[3] <- -log10(temp$value)

names(gtex.counts) <- c('Sample', 'Set.size', 'p.Value')
########################################################################
### Add plot data to counts.
########################################################################

counts.with.bands <- test.t[1,]
counts.with.bands <- rbind(counts.with.bands,rep(9,length(test.t)))
for(row in 2:length(test.t[,1])) {
  counts.with.bands <- rbind(counts.with.bands,test.t[row,])
  counts.with.bands <- rbind(counts.with.bands,rep(9,length(test.t)))
}

test <- gtex.results[1:5,]
test.t <- t(test)







##########################################################
### Removing the empties
##########################################################

sample.pos <- 8

ggplot(data = gtex.counts[((sample.pos-1)*4+1):(((sample.pos-1)*4)+4),], aes(x = Sample, y = Set.size, fill = p.Value)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar()


((sample.pos-1)*4+1):(((sample.pos-1)*4)+4)

### Generate plot from melted data
gtex.plot <- generate.gtex.fig.from.data(gtex.plot.data.melt)

### Display the plot
gtex.plot

### Save plot
tiff(filename=paste0("gtex_",geneset,"_",length,"_",strftime(Sys.time(),"%a%b%d%H%M"),".tiff"), 
     type="cairo",
     units="in", 
     width=14, 
     height=14, 
     pointsize=12,
     res=100)
gtex.plot
dev.off()

### Generate multiple plots in list
start <- 600
end <- 2600
step <- 200

gtex.plots <- list()
for(length in seq(start,end,step)) {
  print(paste("Generating",length(seq(start,end,step)),"plots using",geneset))
  gtex.plots[[length(gtex.plots)+1]] <- generate.gtex.figure(ids.of.interest.full,length = length,FALSE)
  print(paste('plot of',length,'genes complete'))
}

### View plots in list
gtex.plot <- gtex.plots[[8]]
gtex.plot

### Save plots in list
title.split <- strsplit(gtex.plot$labels$title," ")[[12]]
geneset <- title.split[8]
length <- title.split[5]
tiff(filename=paste0("gtex_",geneset,"_",length,"_",strftime(Sys.time(),"%a%b%d%H%M"),".tiff"), 
     type="cairo",
     units="in", 
     width=14, 
     height=14, 
     pointsize=12,
     res=100)
gtex.plot
dev.off()

########################################################################
### Generate tsea data and figures
########################################################################

length = 2222

### Generate the plot
tsea.plot <- generate.tsea.figure(genes.of.interest.full,length,FALSE)

### Generate the gene list from brain
tsea.brain.genes <- 

### Display the plot
tsea.plot

### Save plot
tiff(filename=paste0("tsea_",geneset,"_",length,"_",strftime(Sys.time(),"%a%b%d%H%M"),".tiff"), 
     type="cairo",
     units="in", 
     width=14, 
     height=14, 
     pointsize=12,
     res=100)
tsea.plot
dev.off()

### Define an empty list
brain.genes <- list()

### Generate lists of overlapping genes for the given gene depth
tsea.genes <- candidate.overlap(pSIs = psi.tsea, candidate.genes = genes.of.interest.full[1:length])

### Add brain genes to list
brain.genes[[length(brain.genes)+1]] <- tsea.genes$pSi_0.05$Brain_0.05 ; names(brain.genes)[length(brain.genes)] <- paste0("brain_0.05-",length)
brain.genes[[length(brain.genes)+1]] <- tsea.genes$pSi_0.01$Brain_0.01 ; names(brain.genes)[length(brain.genes)] <- paste0("brain_0.01-",length)
brain.genes[[length(brain.genes)+1]] <- tsea.genes$pSi_0.001$Brain_0.001 ; names(brain.genes)[length(brain.genes)] <- paste0("brain_0.001-",length)



########################################################################
### Calculate overlap p-values
########################################################################

tsea.results <- fisher.iteration(pSIs = psi.tsea, candidate.genes = genes.of.interest)
csea.results <- fisher.iteration(pSIs = psi.csea, candidate.genes = genes.of.interest)
gtex.results <- fisher.iteration(pSIs = psi.gtex, candidate.genes = ids.of.interest)

########################################################################
### Report significant genes
########################################################################

tsea.genes <- candidate.overlap(pSIs = psi.tsea, candidate.genes = genes.of.interest)

brain.genes <- list()

brain.genes[[length(brain.genes)+1]] <- tsea.genes$pSi_0.05$Brain_0.05
brain.genes[[length(brain.genes)+1]] <- tsea.genes$pSi_0.01$Brain_0.01
brain.genes[[length(brain.genes)+1]] <- tsea.genes$pSi_0.001$Brain_0.001


########################################################################
### Save p-value results
########################################################################

timestamp <- strftime(Sys.time(),"%a%b%d%H%M")
setwd("z://Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/pSI_results/")

write.csv(tsea.results, paste0("tsea_",geneset,"_",length,"_",timestamp,".csv"))
write.csv(csea.results, paste0("csea_",geneset,"_",length,"_",timestamp,".csv"))
write.csv(gtex.results, paste0("gtex_",geneset,"_",length,"_",timestamp,".csv"))

########################################################################
### Load p-value results
########################################################################

#setwd("z://Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/pSI_results/")

#tsea.results <- read.csv("tsea_")
#csea.results <- read.csv("csea_")
#gtex.results <- read.csv("gtex_")

########################################################################
### Visualize results
########################################################################
### Plot bar graphs 

### Format TSEA data
tsea.plotting <- tsea.results[tsea.results$`0.05 - adjusted` != 1,]
names(tsea.plotting) <- c('p0.05','p0.01','p0.001','p1e-4')
tsea.plotting <- -log10(tsea.plotting)
tsea.plotting$tissue <- row.names(tsea.plotting)
tsea.plotting.melt <- melt(tsea.plotting,id.var = "tissue")

### Plot TSEA data
tsea.plot <- ggplot(data = tsea.plotting.melt,aes(x = reorder(tissue, -value), y = value, fill = variable)) +
  geom_bar(position = "dodge",stat="identity") + 
  labs(title = paste("p-value of the overlap of the top",length,"genes in",geneset),
       x = "Tissue",
       y = "Negative log10 p-value") +
  theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer("p-value",palette = 'OrRd', direction = -1)

tsea.plot

### Format CSEA data
csea.plotting <- csea.results[csea.results$`0.05 - adjusted` != 1,]
names(csea.plotting) <- c('p0.05','p0.01','p0.001','p1e-4')
csea.plotting <- -log10(csea.plotting)
csea.plotting$tissue <- row.names(csea.plotting)
csea.plotting.melt <- melt(csea.plotting,id.var = "tissue")

### Plot CSEA data
csea.plot <- ggplot(data = csea.plotting.melt,aes(x = reorder(tissue, -value), y = value, fill = variable)) +
  geom_bar(position = "dodge",stat="identity") + 
  labs(title = paste("p-value of the overlap of the top",length,"genes in",geneset),
       x = "Tissue",
       y = "Negative log10 p-value") +
  theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer("p-value",palette = 'OrRd', direction = -1)

csea.plot

### Format gtex data
gtex.plotting <- gtex.results[gtex.results$`0.01 - adjusted` != 1,]
names(gtex.plotting) <- c('p0.05','p0.01','p0.001','p1e-4')
gtex.plotting <- gtex.plotting[c('p0.01','p0.001','p1e-4')]
gtex.plotting <- -log10(gtex.plotting)
gtex.plotting$tissue <- row.names(gtex.plotting)
gtex.plotting[order(gtex.plotting$p0.05,decreasing = TRUE),]
gtex.plotting.melt <- melt(gtex.plotting,id.var = "tissue")

### Plot gtex data
gtex.plot <- ggplot(data = gtex.plotting.melt,aes(x = reorder(tissue, -value), y = value, fill = variable)) +
  geom_bar(position = "dodge",stat="identity") + 
  labs(title = paste("Overlap of the top",length,"genes in",geneset,"and gtex"),
       x = "Tissue",
       y = "Negative log10 p-value") +
  theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer("p-value",palette = 'Reds')

gtex.plot

########################################################################
### Output plots
########################################################################

setwd("z://Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/pSI_results/")

tiff(filename=paste0("tsea_",geneset,"_",length,"_",strftime(Sys.time(),"%a%b%d%H%M"),".tiff"), 
     type="cairo",
     units="in", 
     width=14, 
     height=14, 
     pointsize=12, 
     res=100)
tsea.plot
dev.off()

tiff(filename=paste0("csea_",geneset,"_",length,"_",strftime(Sys.time(),"%a%b%d%H%M"),".tiff"), 
     type="cairo",
     units="in", 
     width=14, 
     height=14, 
     pointsize=12, 
     res=100)
csea.plot
dev.off()







###################################################
### Generate overlap lists
###################################################

test <- candidate.overlap(pSIs=psi.tsea,candidate.genes=genes.of.interest)

test <- candidate.overlap(brain_genes_0.05,genes.of.interest)



