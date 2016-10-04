### Binding site identification
### Andrew R Gross, 2016-09-26

########################################################################
### Header
########################################################################

library(MotifDb)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)
library(org.Sc.sgd.db)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

########################################################################
### Get Pos. Freq. Matrices
########################################################################

### Report the number of detected binding sequences
(RELA <- query(MotifDb,"RELA"))
(TP53 <- query(MotifDb,"TP53"))

### Assign binding sequences to variables
RELA.both <- RELA[[1]]
TP53.core <- TP53[[1]]
TP53.14 <- TP53[[2]]

### Visualize binding site plots
seqLogo(RELA.both)  # Does this look phosphoralated??
seqLogo(TP53.core)
seqLogo(TP53.14)

### Format pfm as integers
pfm.rela <- round(100*RELA.both)
pfm.tp53.core <- round(100*TP53.core)
pfm.tp53.14 <- round(100*TP53.14)

pfm.list <- list(pfm.rela,pfm.tp53.core,pfm.tp53.14)

########################################################################
### Get Promoter regions
########################################################################

### Specify genes of interest
genes <- c("SCO2","TFAM","POLRMT","CYB5A")
#genes.other <- c("DLAT","MDC1","IL10","QARS","ZXDA")
#genes <- c("DAL1", "DAL2", "DAL4", "DAL5", "DAL7", "DAL80", "GAP1")

### Declare pos. controls
genes <- c('IL1A','IL1B','TNF','IL6')  # RELA positive controls   ,'IL-8','MCP-1' not found
genes <- c('CDKN1A','GADD45A','GADD45B','GADD45G','PERP','BAX')   # TP53 Positive controls

### Generate a list of GRange objects from the target organism
x <- org.Hs.egSYMBOL2EG         ## Objects in this package can be accessed using the select() interface from the AnnotationDbi package. See ?select for details.
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

### Format entrez gene IDs to call GRange objects from full list
gene.positions <- xx[genes]
gene.positions <- as.vector(unlist(gene.positions))

### Select the GRange objects corresponding to genes of interest
transcript.seqs <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene") [gene.positions]
promoter.size <- 2000
promoter.seqs <- getPromoterSeq(transcript.seqs,Hsapiens, upstream=promoter.size, downstream=0)
promoter.seqs <- unlist(promoter.seqs)
### Select one example of each unique promoter
promoter.names.unique <- unique(names(promoter.seqs))
unique.promoters <- match(promoter.names.unique,names(promoter.seqs))
promoter.seqs <- promoter.seqs[unique.promoters]

########################################################################
### Get Random Promoter regions
########################################################################

rand.positions <- xx[sample.int(length(xx),size = 10)]
rand.positions <- as.vector(unlist(rand.positions))

samples <- 50
all.transcript.seqs <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene")
rand.transcript.seqs <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene") [sample.int(23459,size = samples)]
rand.promoter.seqs <- getPromoterSeq(rand.transcript.seqs,Hsapiens, upstream=promoter.size, downstream=0)
rand.promoter.seqs <- unlist(rand.promoter.seqs)

### Select for unique values
rand.seq.uniq.name <- unique(names(rand.promoter.seqs))
rand.promoter.seqs <- rand.promoter.seqs[match(rand.seq.uniq.name,names(rand.promoter.seqs))]

########################################################################
### Identify pfm matches in promoter sequences
########################################################################

pfm.rela
pfm.tp53.core
pfm.tp53.14
pcm.dal80.jaspar

matchPWM(pfm.tp53.14, promoter.seqs[[2]], "70%") ->test
matchPWM(pfm.rela, promoter.seqs[[2]], "70%")
matchPWM(pfm.list[[1]], promoter.seqs[[2]], "70%")

pwm.hits <- sapply(promoter.seqs, function(pseq) matchPWM(pfm.list[[1]], pseq, min.score="70%"))
test <- sapply(pwm.hits, length)

########################################################################
### Compare all transcription factors and targets
########################################################################
column.pos <- 1
sapply(promoter.seqs, function(pseq) matchPWM(pfm.list[[column.pos]], pseq, min.score="60%"))
#
### Calculate columns of hits
#column.list <- list()
gene.number <- length(pfm.list)
binding.hits <- data.frame(gene = names(promoter.seqs))

for(column.pos in 1:gene.number) {
  temp.list <- sapply(promoter.seqs, function(pseq) matchPWM(pfm.list[[column.pos]], pseq, min.score="60%"))
  temp.vector <- sapply(temp.list,  length)
  binding.hits[ncol(binding.hits)+1] <- temp.vector
  #column.list[[column.pos]] <- temp.vector
}
names(binding.hits) <- c("genes","RELA","TP53.co","TP53.14")


########################################################################
### Compare TF binding to random targets
########################################################################
column.pos <- 1
rand.results <- sapply(rand.promoter.seqs, function(pseq) matchPWM(pfm.list[[column.pos]], pseq, min.score="85%"))
rand.results <- sapply(promoter.seqs, function(pseq) matchPWM(pfm.list[[column.pos]], pseq, min.score="60%"))

hit.numbers <- c()
closest.prox <- c()

for(item in 1:length(rand.results)){
  # Calculate hit numbers
  hit.numbers <- c(hit.numbers,length(rand.results[[item]]))
  # Calculate most proximal binding location
  
  new.pos <- start(rand.results[[item]])[1]
  if(is.na(new.pos)) {new.pos <- 9999}
  #print(new.pos)
  closest.prox <- c(closest.prox,new.pos)
}

### Order results
hit.numbers <- sort(hit.numbers,decreasing = TRUE)
closest.prox <- sort(closest.prox)

### Add to dataframe
#rand.hits.mat <- matrix(hit.numbers, nrow = samples)
#rand.pos.mat <- matrix(closest.prox, nrow = samples)

rand.hits.mat <- cbind(rand.hits.mat,hit.numbers)
rand.pos.mat <- cbind(rand.pos.mat,closest.prox)

#names(rand.hits.mat) <- c("60%","70%","80%","85%")
#names(rand.pos.mat) <- c("60%","70%","80%","85%")

########################################################################
### Compare TF binding to selected genes
########################################################################
column.pos <- 3
results <- sapply(promoter.seqs, function(pseq) matchPWM(pfm.list[[column.pos]], pseq, min.score="85%"))

hit.numbers <- c()
closest.prox <- c()

for(item in 1:length(results)){
  # Calculate hit numbers
  hit.numbers <- c(hit.numbers,length(results[[item]]))
  # Calculate most proximal binding location
  
  new.pos <- start(results[[item]])[1]
  if(is.na(new.pos)) {new.pos <- 9999}
  #print(new.pos)
  closest.prox <- c(closest.prox,new.pos)
}

### Order results
hit.numbers <- sort(hit.numbers,decreasing = TRUE)
closest.prox <- sort(closest.prox)

### Add to dataframe
#rand.hits.mat <- matrix(hit.numbers, nrow = length(promoter.seqs))
#rand.pos.mat <- matrix(closest.prox, nrow = length(promoter.seqs))

rand.hits.mat <- cbind(rand.hits.mat,hit.numbers)
rand.pos.mat <- cbind(rand.pos.mat,closest.prox)

#names(rand.hits.mat) <- c("60%","70%","80%","85%")
#names(rand.pos.mat) <- c("60%","70%","80%","85%")



## ------------------------------------------------------------------------
grl <- transcriptsBy(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, by="gene") [orfs]

## ------------------------------------------------------------------------
promoter.seqs <- getPromoterSeq(grl, Scerevisiae, upstream=1000,downstream=0)



library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
e2f3 <- "1871" # entrez geneID for a cell cycle control transcription
# factor, chr6 on the plus strand
transcriptCoordsByGene.GRangesList.full <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene")
transcriptCoordsByGene.GRangesList <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene") [e2f3]
transcriptCoordsByGene.GRangesList <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene") ["7019"]

# a GrangesList of length one, describing three transcripts
promoter.seqs <- getPromoterSeq (transcriptCoordsByGene.GRangesList,
                                 Hsapiens, upstream=10, downstream=0)





########################################################################
### Identify matches
########################################################################

query(MotifDb, "DAL80")   
pfm.dal80.jaspar <- query(MotifDb,"DAL80")[[1]]
seqLogo(pfm.dal80.jaspar)
dal1 <- "YIR027C"
chromosomal.loc <- transcriptsBy(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, by="gene") [dal1]
promoter.dal1 <- getPromoterSeq(chromosomal.loc, Scerevisiae, upstream=1000, downstream=0)
pcm.dal80.jaspar <- round(100 * pfm.dal80.jaspar)
matchPWM(pcm.dal80.jaspar, unlist(promoter.dal1)[[1]], "90%")

## ------------------------------------------------------------------------
query(MotifDb,"DAL80")

#NFKB1 <- query(MotifDb,"NFKB1")


## ------------------------------------------------------------------------
dal80.jaspar <- query(MotifDb,"DAL80")[[1]]
dal80.scertf <-query(MotifDb,"DAL80")[[2]]
seqLogo(dal80.jaspar)
seqLogo(dal80.scertf)

pfm.dal80.jaspar <- new("pfm", mat=query(MotifDb, "dal80")[[1]], name="DAL80-JASPAR")
pfm.dal80.scertf <- new("pfm", mat=query(MotifDb, "dal80")[[2]], name="DAL80-ScerTF")
plotMotifLogoStack(DNAmotifAlignment(c(pfm.dal80.scertf, pfm.dal80.jaspar)))

## ---- dev="jpeg"---------------------------------------------------------
pfm.dal80.jaspar <- new("pfm", mat=query(MotifDb, "dal80")[[1]], 
                        name="DAL80-JASPAR")
pfm.dal80.scertf <- new("pfm", mat=query(MotifDb, "dal80")[[2]], 
                        name="DAL80-ScerTF")
plotMotifLogoStack(DNAmotifAlignment(c(pfm.dal80.scertf, pfm.dal80.jaspar)))

## ------------------------------------------------------------------------
query(MotifDb, "gat1")

## ---- dev="jpeg"---------------------------------------------------------
pfm.gat1.jaspar = new("pfm", mat=query(MotifDb, "gat1")[[1]], 
                      name="GAT1-JASPAR")
pfm.gat1.scertf = new("pfm", mat=query(MotifDb, "gat1")[[2]], 
                      name="GAT1-ScerTF")
pfm.gat1.uniprobe = new("pfm", mat=query(MotifDb, "gat1")[[3]], 
                        name="GAT1-UniPROBE")
plotMotifLogoStack(c(pfm.gat1.uniprobe, pfm.gat1.scertf, pfm.gat1.jaspar))

## ------------------------------------------------------------------------
pfm.dal80.scertf <- query(MotifDb, "dal80")[[2]]
pcm.dal80.scertf <- round(100 * pfm.dal80.scertf)

pfm.gat1.jaspar <- query(MotifDb, "gat1")[[1]]
pcm.gat1.jaspar <- round(100 * pfm.gat1.jaspar)

pfm.gat1.scertf <- query(MotifDb, "gat1")[[2]]
pcm.gat1.scertf <- round(100 * pfm.gat1.scertf)

## ------------------------------------------------------------------------
genes <- c("DAL1", "DAL2", "DAL4", "DAL5", "DAL7", "DAL80", "GAP1")
orfs <- as.character(mget(genes, org.Sc.sgdCOMMON2ORF))

## ------------------------------------------------------------------------
grl <- transcriptsBy(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, by="gene") [orfs]
grl2 <- transcriptsBy(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, by="gene")

## ------------------------------------------------------------------------
promoter.seqs <- getPromoterSeq(grl, Scerevisiae, upstream=1000,
                                downstream=0)

## ------------------------------------------------------------------------
pfm.dal80.scertf

## ------------------------------------------------------------------------
print (class(promoter.seqs))
promoter.seqs <- unlist(promoter.seqs)
print (class(promoter.seqs))

matchPWM(pcm.dal80.scertf, promoter.seqs[[1]], "90%")

## ------------------------------------------------------------------------
pwm.hits <- sapply(promoter.seqs, function(pseq) matchPWM(pcm.dal80.scertf, pseq, min.score="90%"))

## ------------------------------------------------------------------------
dal80.scertf.hits <- sapply(promoter.seqs, function(pseq) matchPWM(pcm.dal80.scertf, pseq, min.score="90%"))
gat1.scertf.hits  <- sapply(promoter.seqs, function(pseq) matchPWM(pcm.gat1.scertf, pseq, min.score="90%"))
gat1.jaspar.hits  <- sapply(promoter.seqs, function(pseq) matchPWM(pcm.gat1.jaspar, pseq, min.score="90%"))

## ------------------------------------------------------------------------
dal80.scertf <- sapply(dal80.scertf.hits, length)
gat1.jaspar  <- sapply(gat1.jaspar.hits,  length)
gat1.scertf  <- sapply(gat1.scertf.hits,  length)

## ------------------------------------------------------------------------
tbl.gata     <- data.frame(gene=genes, dal80.scertf, gat1.jaspar, gat1.scertf)



