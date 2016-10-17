### Tissues Specific Expression Analysis, Andrew R Gross, 2016-08-26
### Using a package developed by the Dougherty lab, calculate the uniqueness of a set of genes and the probability
### of another set containing those genes.

########################################################################
### Header
########################################################################
library(pSI)
load("C:/Users/grossar/Bioinform/DATA/rna_seq/Dougherty/pSI.data/data/human.rda")
help(pSI)

ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

########################################################################
### Functions
########################################################################

########################################################################
### Data input
########################################################################

### Human Brain - Microarray
setwd('C://Users/grossar/Bioinform/DATA/Allen Brain Atlas/Human Brain/Normalized microarray datasets/')
human.micro.61 <- read.csv('normalized_microarray_donor9861/MicroarrayExpression.csv')
samples.61 <- read.csv('normalized_microarray_donor9861/SampleAnnot.csv')
ontology.61 <- read.csv('normalized_microarray_donor9861/Ontology.csv')

human.micro.21 <- read.csv('normalized_microarray_donor10021/MicroarrayExpression.csv')
samples.61 <- read.csv('normalized_microarray_donor10021/SampleAnnot.csv')
ontology.61 <- read.csv('normalized_microarray_donor10021/Ontology.csv')

human.micro.76 <- read.csv('normalized_microarray_donor12876/MicroarrayExpression.csv')
human.micro.80 <- read.csv('normalized_microarray_donor14380/MicroarrayExpression.csv')
human.micro.96 <- read.csv('normalized_microarray_donor15496/MicroarrayExpression.csv')
human.micro.97 <- read.csv('normalized_microarray_donor15697/MicroarrayExpression.csv')

### Developing Brain - Microarray
setwd('C://Users/grossar/Bioinform/DATA/Allen Brain Atlas/Developing Brain/')
developing.micro.66 <- read.csv('lmd_matrix_12566/expression_matrix.csv', header = FALSE, row.names = 1)
developing.micro.90 <- read.csv('lmd_matrix_12690/expression_matrix.csv', header = FALSE, row.names = 1)
developing.micro.40 <- read.csv('lmd_matrix_12840/expression_matrix.csv', header = FALSE, row.names = 1)
developing.micro.51 <- read.csv('lmd_matrix_14751/expression_matrix.csv', header = FALSE, row.names = 1)

dev.genes.66 <- read.csv('lmd_matrix_12566/rows_metadata.csv')

dev.micro.66.unq <- avg.redundant.probes(developing.micro.66, dev.genes.66$entrez_id)



genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name','entrezgene'), filters='ensembl_gene_id', values=ids.of.interest.full, mart=ensembl)


length(unique(dev.genes.66$entrez_id))

avg.redundant.probes <- function(dataframe, ids) {   # Average redundant rows into a new dataframe
  dataframe <- dataframe[-]
  output.dataframe <- dataframe[1,]   # Declare new dataframe with first row copied from input dataframe
  number.of.probes <- length(unique(ids))   # Declare the number of unique entrez ids
  previous.id <- 0   # Declare the 'previous id'
  for(row.position in 2:length(ids)) {   # begin looping through input dataframe
    current.id <- ids[row.position]   # Declare the ID of the current row
    if(current.id > previous.id) {     # If the current id is different from the previous one...
      output.dataframe[nrow(output.dataframe)+1,] <- dataframe[row.position,]   # ...copy it to the new dataframe.
    }
    if(current.id == previous.id) {   # If the current id is the same as the previous one...
      rows.to.average <- rbind(output.dataframe[nrow(output.dataframe),],dataframe[row.position,])   # Generate a dataframe containing the new and old rows
      averaged.row <- apply(rows.to.average,2,mean)  # And make a new row from the average of the two
      output.dataframe[nrow(output.dataframe),] <- averaged.row  # replace the last row with the new averaged row.
    }
    previous.id <- current.id
  }
  row.names(output.dataframe) <- unique(ids)
  return(output.dataframe[1:ncol(output.dataframe)-1])
}

gene.name.temp <- dev.genes.66$gene_name[nas]
temp.expression.data <- developing.micro.66[nas,]


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
gtex.counts <- read.csv('Z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_gtex_counts_10_03_16.csv', row.names = 1)

########################################################################
### pSI determination
########################################################################

#pSI.input <- sample.data$pSI.input
#pSI.input <- references[2:ncol(references)]
#pSI.output <- specificity.index(pSI.input)
#pSI.thresholds <- pSI.list(pSI.output)
#psi.counts <- pSI.count(psi.gtex)

########################################################################
### Subsetting
########################################################################

### Subset metadata files -- Metadata files contain data used for subsetting data sets and formatting plots
metadata.sam.HT <- na.omit(metadata.sam[metadata.sam$Type == "HT",])       # all HT samples (iPSC + adult)
metadata.sam.iHT <- metadata.sam[metadata.sam$Group == 'iHT',]             # iPSC HT samples
metadata.sam.aHT <- metadata.sam[metadata.sam$Source == "Adult",]          # Adult HT samples
metadata.sam.iPSC <- metadata.sam[metadata.sam$Source == "iPSC",]          # all iPSC samples (iHT + iMN)
metadata.sam.MN <- na.omit(metadata.sam[metadata.sam$Type == "MN",])       # Motor neuron samples
metadata.sam.HT.f <- metadata.sam.HT[metadata.sam.HT$Sex =='F',]           # female HT samples
metadata.sam.iHT.f <- metadata.sam.iHT[metadata.sam.iHT$Sex =='F',]        # female iHT samples

### Subset expression data
HT.data <- TPMdata[match(row.names(metadata.sam.HT),names(TPMdata))]       # All hypothalamus samples (iPSC + adult)
iHT.data <- TPMdata[match(row.names(metadata.sam.iHT),names(TPMdata))]     # iPSC hypothalamus samples
aHT.data <- TPMdata[match(row.names(metadata.sam.aHT),names(TPMdata))]     # Adult hypothalamus samples
iPSC.data <- TPMdata[match(row.names(metadata.sam.iPSC),names(TPMdata))]   # all iPSC samples (iHT + iMN)
MN.data <- TPMdata[match(row.names(metadata.sam.MN),names(TPMdata))]       # Motor Neuron samples
iHT.f.data <- TPMdata[match(row.names(metadata.sam.iHT.f),names(TPMdata))] # female HT samples
HT.f.data <- TPMdata[match(row.names(metadata.sam.HT.f),names(TPMdata))]   # female iHT samples

### Join all data sets into a list from which to select
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
geneset <- 2  # Specifies the data set to use
genes.of.interest.full <- genelist[[geneset]][,1] # Declares the genes to use
ids.of.interest.full <- row.names(genelist[[geneset]]) # Declares the IDs of the genes to use
geneset <- names(genelist[geneset]) # Names the geneset for use in plot titles
print(geneset) # Reports the geneset in use

########################################################################
### Generate gtex plot data
########################################################################

length <- 1200

ids.of.interest <- ids.of.interest.full[1:length] # Subselects based on the specified depth
### Calculate p-values ############################################
gtex.results <- fisher.iteration(pSIs = psi.gtex, candidate.genes = ids.of.interest)

### Format gtex data ##############################################
#absent.samples <- row.names(gtex.results[gtex.results$`0.05 - adjusted` == 1,])
names(gtex.results) <- c('p0.05','p0.01','p0.001','p1e-4')

########################################################################
########################################################################
### Export/Import p-value table
########################################################################
########################################################################

setwd("z://Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/pSI_results/")
write.csv(gtex.results, paste0("gtex_",geneset,"_",length,"_",strftime(Sys.time(),"%a%b%d%H%M"),".csv")); print(paste('Wrote',paste0("gtex_",geneset,"_",length,"_",strftime(Sys.time(),"%a%b%d%H%M"))))
gtex.results <- read.csv('gtex_iHT_800_WedOct051352.csv',row.names = 1)

########################################################################
### Format counts and results
########################################################################

### Remove empty rows
#gtex.results.trim <- gtex.results
gtex.results.trim <- gtex.results[gtex.results$p0.05 < 0.1,]
gtex.counts.trim <- gtex.counts[row.names(gtex.results.trim),]
gtex.counts.trim$p1e.4 <- gtex.counts.trim$p1e.4+1

### Reverse row order
gtex.counts.trim <- gtex.counts.trim[c(4,3,2,1)]
gtex.results.trim <- gtex.results.trim[c(4,3,2,1)]

### Melt and join
gtex.plotting <- melt(t(gtex.counts.trim))
gtex.plotting[4] <- melt(t(gtex.results.trim))[3]
gtex.plotting[5] <- gtex.plotting[3]
gtex.plotting <- gtex.plotting[c(2,3,4,5)]
names(gtex.plotting) <- c('Sample', 'Set.size', 'p.Value', 'bands')

### Rescale and rename
gtex.plotting[2] <- log10(gtex.plotting[2]+1)
gtex.plotting[3] <- -log10(gtex.plotting[3])
#gtex.plotting[4] <- log10(gtex.plotting[4]+1)

### Add bands between bands
gtex.plotting <- add.bands.between(gtex.plotting,0.25)

### Calculate maximum bullseye size for plot limits
plot.max.lim <- ceiling(max(rowSums(log10(gtex.counts.trim))))

########################################################################
### Plot one bullseye
########################################################################

### Declare the gene of interest by giving its position
sample.pos <- 2

### Declare a subsetted dataframe using the specified gene of interest's position
gtex.plotting.sub <- gtex.plotting[(((sample.pos-1)*8)+1):(((sample.pos-1)*8)+8),]

### Generate plot
ggplot(data = gtex.plotting.sub, aes(x = Sample, y = Set.size, fill = p.Value)) +
  geom_bar(width = 1, stat = "identity") + 
  scale_y_continuous(limits = c(0,12)) +
  xlab(gtex.plotting.sub[1][1,]) +
  theme(axis.title.x = element_text(face="bold", size=14,margin =margin(0,0,10,0)),
        axis.title.y = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"), axis.text = element_text(size = 12),
        panel.background = element_rect(fill = 'white')) +  coord_polar() +
  scale_fill_distiller(palette = 'YlOrRd', direction = 1, limits=c(1, 50), oob = censor, na.value = 'gray95') +
  geom_bar(aes(y = bands), width = 1, stat = 'identity', fill = c(alpha('black',0),'black',alpha("black",0),'black',alpha("black",0),'black',alpha("black",0),'black'))
  
########################################################################
### Plot tissue types in batch for faceting
########################################################################

plot.list <- generate.bullseyes(gtex.plotting, plot.max.lim)

########################################################################
### Generate faceted plots
########################################################################
print(paste('There are', length(plot.list), 'tissue types available for plotting'))

### Specify bullseyes to plot
p.l <- plot.list

n <- length(p.l); nCol <- floor(sqrt(n))

(tiled.figure <- do.call("grid.arrange", c(p.l, ncol=nCol)))




(tiled.figure <- grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],
                              ncol = 4, top = title))


##########################################################
### Saving faceted plots
##########################################################

title.split <- strsplit(gtex.plot$labels$title," ")[[12]]
geneset <- title.split[8]
length <- title.split[5]

getwd()
png(filename=paste0("gtex_",geneset,"_",length,"_",strftime(Sys.time(),"%a%b%d%H%M"),".png"), 
     type="cairo",
     units="in", 
     width=14, 
     height=14, 
     pointsize=12,
     res=300)
do.call("grid.arrange", c(p.l, ncol=nCol))
dev.off()



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


########################################################################
### Scratch work
########################################################################


(tiled.figure<-grid.arrange(p.l[[1]],p.l[[2]],p.l[[3]],p.l[[4]],p.l[[5]],
                            p.l[[6]],p.l[[7]],p.l[[8]],p.l[[9]],p.l[[10]],
                            p.l[[11]],p.l[[12]],p.l[[13]],p.l[[14]],p.l[[15]],
                            ncol = 5, top = title))
