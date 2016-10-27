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
(RELA <- query(MotifDb,"RELA"))  # NFkB
(TP53 <- query(MotifDb,"TP53"))  # p53
(NFKB <- query(MotifDb,"NFKB"))  #

### Assign binding sequences to variables
RELA.both <- RELA[[1]]
TP53.core <- TP53[[1]]
TP53.14 <- TP53[[2]]
NFKB <- NFKB[[4]]

### Visualize binding site plots
seqLogo(RELA.both)  # Does this look phosphoralated??
seqLogo(TP53.core)
seqLogo(TP53.14)
seqLogo(NFKB)

### Format pfm as integers
pfm.rela <- round(100*RELA.both)
pfm.tp53.core <- round(100*TP53.core)
pfm.tp53.14 <- round(100*TP53.14)
pfm.nfkb <- round(100*NFKB)

### Create list containing all pfm
pfm.list <- list(pfm.rela,pfm.tp53.core,pfm.tp53.14,pfm.nfkb)
names(pfm.list) <- c("RELA","TP53.core","TP53.14",'NFKB')

########################################################################
### Get Promoter regions
########################################################################

### Specify genes of interest
genes.unknown1 <- c("SCO2","TFAM","POLRMT","CYB5A")
genes.unknown2 <- c("CHUK")
genes.unknown <- c(genes.unknown1, genes.unknown2)

### Declare pos. controls
genes.pos1 <- c('IL1A','IL1B','TNF','IL6')  # RELA positive controls   ,'IL-8','MCP-1' not found
genes.pos2 <- c('CDKN1A','GADD45A','GADD45B','GADD45G','PERP','BAX')   # TP53 Positive controls
genes.pos <- c(genes.pos1,genes.pos2)

### Declare negative controls
genes.neg <- c()
#genes.other <- c("DLAT","MDC1","IL10","QARS","ZXDA")
#genes <- c("DAL1", "DAL2", "DAL4", "DAL5", "DAL7", "DAL80", "GAP1")

### Create full list of genes to query
genes <- c(genes.unknown, genes.pos, genes.neg)
gene.type <- c(rep('unknown',length(genes.unknown)),rep('pos.',length(genes.pos)),rep('neg',length(genes.neg)))

### Generate a list of GRange objects from the target organism
x <- org.Hs.egSYMBOL2EG         # Objects in this package can be accessed using the select() interface from the AnnotationDbi package. See ?select for details.
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

### Generate a list of Entrez IDs from the list of specified genes
gene.positions <- as.vector(unlist(xx[genes]))

### Select the GRange objects corresponding to genes of interest
transcript.seqs <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene") [gene.positions] # Generate GRange objects
promoter.size <- 2000  # Specify promoter size
promoter.seqs <- getPromoterSeq(transcript.seqs,Hsapiens, upstream=promoter.size, downstream=0)  # Call promoter with the getPromoterSeq() comand
promoter.seqs <- unlist(promoter.seqs)

########################################################################
### Get Random Promoter regions
########################################################################

### Call protomter sequences for random transcripts
samples <- 50
rand.transcript.seqs <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene") [sample.int(23459,size = samples)]
rand.promoter.seqs <- unlist(getPromoterSeq(rand.transcript.seqs,Hsapiens, upstream=promoter.size, downstream=0))

### Select for unique values
rand.seq.uniq.name <- unique(names(rand.promoter.seqs))
rand.promoter.seqs <- rand.promoter.seqs[match(rand.seq.uniq.name,names(rand.promoter.seqs))]


########################################################################
### Generate dataframe of hit number results
########################################################################

### Define columns for dataframe
trans.factor <- c()
for(tf.pos in 1:length(pfm.list)){trans.factor <- c(trans.factor, rep(names(pfm.list)[tf.pos],length(genes)))}
gene <- rep(genes,length(pfm.list))
type <- rep(gene.type,length(pfm.list))

### Generate a dataframe to populate
hit.numbers.df <- data.frame(trans.factor, gene, type)

### Define specificity levels to query
specificity.levels <- c('60%', '70%', '80%', '85%', '90%', '95%')

### Populate dataframe with hit numbers
for(spec.level in specificity.levels){
  print(spec.level)
  hit.numbers.for.spec.lvl <- c()
  ### Calculate the number of hits for each gene against each binding site
  for(pfm in pfm.list){
    results <- sapply(promoter.seqs, function(pseq) matchPWM(pfm, pseq, min.score = spec.level)) # Generate list of matches in unfamiliar format
    unique.gene.positions <- unique(names(results))
    ### For the given motif and specificity find the hit numbers of all unique gene names
    hit.numbers.for.tf <- c()
    for(current.gene in unique.gene.positions){
      sub.results <- results[grep(current.gene,names(results))]  # Subset the results into a list containing variants of the same transcript
      sub.hit.numbers <- c()
      for(item in sub.results){      # Calculate the number of hits for each of the transcripts and add them to a vector
        sub.hit.numbers <- c(sub.hit.numbers,length(item))
      }
      hit.number <- max(sub.hit.numbers)        # Calculate the greatest number of hits and
      hit.numbers.for.tf <- c(hit.numbers.for.tf, hit.number) # add that to the vector of hits for other genes
    }
    hit.numbers.for.spec.lvl <- c(hit.numbers.for.spec.lvl, hit.numbers.for.tf)
  }
  hit.numbers.df <- cbind(hit.numbers.df,hit.numbers.for.spec.lvl)
}
names(hit.numbers.df) <- c('trans.factor', 'gene', 'type', 'hits.60%', 'hits.70%', 'hits.80%', 'hits.85%', 'hits.90%', 'hits.95%')

print(hit.numbers.df)

setwd('z:/Uthra/EDC paper/Binding sites/')
#write.csv(nearest.matches, paste0("nearest.matches_",strftime(Sys.time(),"%a%b%d%H%M"),".csv")); print(paste('Wrote',paste0("nearest.matches_",strftime(Sys.time(),"%a%b%d%H%M"))))

########################################################################
### Generate dataframe of nearest position results
########################################################################

### Define columns for dataframe
trans.factor <- c()
for(tf.pos in 1:length(pfm.list)){trans.factor <- c(trans.factor, rep(names(pfm.list)[tf.pos],length(genes)))}
gene <- rep(genes,length(pfm.list))
type <- rep(gene.type,length(pfm.list))

### Generate a dataframe to populate
nearest.matches.df <- data.frame(trans.factor,gene,type)

### Define specificity levels to query
specificity.levels <- c('60%', '70%', '80%', '85%', '90%', '95%')

### Populate dataframe with hit numbers
for(spec.level in specificity.levels){
  print(spec.level)
  nearest.matches.for.spec.lvl <- c()
  ### Calculate the number of hits for each gene against each binding site
  for(pfm in pfm.list){
    results <- sapply(promoter.seqs, function(pseq) matchPWM(pfm, pseq, min.score = spec.level)) # Generate list of matches in unfamiliar format
    unique.gene.positions <- unique(names(results))
    ### For the given motif and specificity find the hit numbers of all unique gene names
    nearest.matches.for.tf <- c()
    for(current.gene in unique.gene.positions){
      sub.results <- results[grep(current.gene,names(results))]  # Subset the results into a list containing variants of the same transcript
      sub.nearest.matches <- c()
      for(item in sub.results){      # Calculate the number of hits for each of the transcripts and add them to a vector
        nearest.match <- start(item)[1]
        if(is.na(nearest.match)) {nearest.match <- '-'}
        sub.nearest.matches <- c(sub.nearest.matches,nearest.match)
    }
      nearest.match <- min(sub.nearest.matches)        # Calculate the greatest number of hits and
      nearest.matches.for.tf <- c(nearest.matches.for.tf, nearest.match) # add that to the vector of hits for other genes
    }
    nearest.matches.for.spec.lvl <- c(nearest.matches.for.spec.lvl, nearest.matches.for.tf)
  }
  nearest.matches.df <- cbind(nearest.matches.df,nearest.matches.for.spec.lvl)
}
names(nearest.matches.df) <- c('trans.factor', 'gene', 'type', 'pos.60%', 'pos.70%', 'pos.80%', 'pos.85%', 'pos.90%', 'pos.95%')

print(nearest.matches.df)

setwd('z:/Uthra/EDC paper/Binding sites/')
#write.csv(nearest.matches, paste0("nearest.matches_",strftime(Sys.time(),"%a%b%d%H%M"),".csv")); print(paste('Wrote',paste0("nearest.matches_",strftime(Sys.time(),"%a%b%d%H%M"))))




########################################################################
### Compare TF binding to random targets
########################################################################
column.pos <- 3
rand.results <- sapply(rand.promoter.seqs, function(pseq) matchPWM(pfm.list[[column.pos]], pseq, min.score="85%"))

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

### Calculate minimum, max, and median values
summary(rand.pos.mat)
summary(rand.hits.mat)











