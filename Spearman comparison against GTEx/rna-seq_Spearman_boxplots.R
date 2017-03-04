### Spearman plotting of HT data -- Andrew R Gross -- 2017-02-24

########################################################################
### Header
library(gplots)
library(ggplot2)
library(Hmisc)
library(reshape2)
library(grid)
library("gridExtra")
library("cowplot")

########################################################################
### Functions
spearman.calc <- function(sample, reference.input) {           # Calculate the spearman correlation between the sample and the references
  spearman.results <- data.frame(rep(0,ncol(reference.input))) # Generate empty results table
  row.names(spearman.results) <- names(reference.input)        # Name empty results table
  names(spearman.results) <- names(sample)
  for(ref.tissue.num in 1:ncol(reference.input)) {             # Loop through each reference 
    ref.tissue.data <- reference.input[ref.tissue.num]         # Call the current tissue from the references
    tissue <- names(ref.tissue.data)
    ref.tissue.data.2 <- ref.tissue.data[which(ref.tissue.data[,1]>0),,drop=FALSE] # Filter out missing values from tissue
    genes.present.in.ref <- row.names(ref.tissue.data.2)         # Declare the genes present in the reference
    genes.missing.in.query <- setdiff(genes.present.in.ref, row.names(sample)) # Declare genes in reference missing from sample
    rows.to.add <- data.frame(rep(0,length(genes.missing.in.query)))  # Generate a zero data frame the size of the missing rows 
    row.names(rows.to.add) <- genes.missing.in.query           # Name rows after missing rows
    names(rows.to.add) <- names(sample)                        # Name column the same as the query 
    sample.2 <- rbind(sample,rows.to.add)                        # Use rbind to make a full data frame containing exactly the genes in the reference
    sample.3 <- sample.2[genes.present.in.ref, , drop = FALSE]     # Reorder sample to match reference
    spearman.input <- cbind(sample.3, ref.tissue.data.2)           # Bind sample and reference
    result <- rcorr(as.matrix(spearman.input), type = 'spearman')[[1]] # Perform spearman calculation
    result <- round(result[2] * 100, 1)                        # Round result
    spearman.results[tissue,] <- result                        # Add to results table
  }
  return(spearman.results)
}
spearman.calc.for.multiple.samples <- function(samples, reference.input) {
  ### Initialize empty results table
  spearman.results <- data.frame(rep(0,ncol(reference.input))) # Generate empty results table
  row.names(spearman.results) <- names(reference.input)        # Name empty results table
  
  ### Calculate Spearman correlations
  for(sample.number in 1:ncol(samples)) {                      # Loop through all 8 samples
    spearman.results[sample.number] <- spearman.calc(samples[sample.number],reference.input)
  }
  ### Reorder results
  row.means <- apply(spearman.results,1,mean)                  # Calculate the average score for a tissue across all 8 samples
  spearman.results <- spearman.results[order(row.means, decreasing = TRUE), , drop = FALSE] # Order the results from highest average tissue to lowest
}
### Plot single boxplot from results table function
plot.tissue.match.boxplot <- function(spearman.results, title) {
  
  ### Reshape data for compatibility with geom_boxplot
  new.col <- ncol(spearman.results)+1 
  spearman.results[new.col] <- row.names(spearman.results)                    # Add the tissues as a new row
  names(spearman.results)[new.col] <- 'ref'                                   # Label the new row
  spearman.t <- t(spearman.results[-new.col])                                 # Transpose the data frame, minus new row
  spearman.melt <- melt(spearman.t)                                     # Melt transposed data frame
  names(spearman.melt) <- c("sample","ref","value")                     # Rename columns of melted data frame
  spearman.melt$ref <- factor(spearman.melt$ref, levels = rev(row.names(spearman.results)), ordered = TRUE) # Assign order
  
  ### Generate data frame for color and label data
  spearman.for.color <- spearman.results                                # Duplicate spearman results
  spearman.for.color[-new.col] <- apply(spearman.for.color[-new.col], 1, mean)    # Reassign all values to mean of group
  names(spearman.for.color)[2] = 'color.val'
  spearman.t <- t(spearman.for.color[-new.col])                               # Transpose again
  spearman.melt.color <- melt(spearman.t)                               # Melt again
  spearman.melt$color <- spearman.melt.color$value                      # Copy repeating values as 'color column to main df
  
  g <- ggplot(data = spearman.melt, aes(x = ref, y = -value, fill = color)) +
    geom_boxplot() +
    coord_flip() +
    scale_fill_gradientn(colors = c('red','orange','yellow','white')) +
    scale_y_continuous(limits = c(-100,-30), position = 'right') +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = 'none', panel.background = element_rect(fill = 'grey97'),
          plot.title = element_text(size = 14, face = 'bold', margin = margin(5,0,5,0), hjust = 0.5)) +
    labs(title = row.names(spearman.results)[1],
         x = '',
         y = title) +
    geom_text(data = spearman.for.color, aes(x = ref, y = -color.val+8, fill = color.val, label = ref), 
              hjust = 0)
  return(g)
}
### Plot multiplot with all five filter levels
multiplot.spearman.results.list <- function(spearman.results.list) {
  boxplot.list <- list()
  filter.levels <- c('Full', 'Loose', 'Neutral', 'Tight', 'Supertight')
  level.num = 1
  for (spearman.results in spearman.results.list) {
    filter.level = filter.levels[level.num]
    level.num <- level.num + 1
    boxplot.list[[length(boxplot.list)+1]] <- plot.tissue.match.boxplot(spearman.results, filter.level)
  }
  print('All figures complete')
  multiboxplot <- plot_grid(boxplot.list[[1]], boxplot.list[[2]], boxplot.list[[3]], boxplot.list[[4]], boxplot.list[[5]], labels=c("F", "L", "N", "T", "ST"), ncol = 5, nrow = 1)
  #multiboxplot <- multiplot(boxplot.list[[1]], boxplot.list[[2]], boxplot.list[[3]], boxplot.list[[4]], boxplot.list[[5]], cols = 5)
  return(multiboxplot)
}
########################################################################
### Data input

### Import sample key -- This table includes the ID, Tissue type, Tissue group, and unique name for all 8555 GTEx samples
sample.key <- read.csv('Z:/Data/Andrew/reference_data/gtex/Sample.key.sorted.csv',header=TRUE)

### Import tables containing all tissue for each level -- These tables each contain the averages of each tissue, filtered for expression and sd
gtex.low.sd.loose <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.loose.csv', row.names = 1)
gtex.low.sd.neutral <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.neutral.csv', row.names = 1)
gtex.low.sd.tight <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.tight.csv', row.names = 1)
gtex.low.sd.supertight <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.supertight.csv', row.names = 1)

### Import full references  -- The average expression of each tissue without filtering
references <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.full.csv',header=TRUE,row.names=1)

### Import known samples  -- A table containing transcriptomes of 8 samples from each tissue in the GTEx collection
reference.transcriptomes <- read.csv('Z:/Data/Andrew/reference_data/gtex/individual-ref-transcriptomes-for-testing.csv', row.names = 1) # ~10 seconds

### iPSC and adult HT data
TPMdata <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)

########################################################################
### Format

### Remove Description column
ref.full <- references[2:ncol(references)]
ref.loose <- gtex.low.sd.loose[2:ncol(gtex.low.sd.loose)]
ref.neutral <- gtex.low.sd.neutral[2:ncol(gtex.low.sd.neutral)]
ref.tight <- gtex.low.sd.tight[2:ncol(gtex.low.sd.tight)]
ref.supertight <- gtex.low.sd.supertight[2:ncol(gtex.low.sd.supertight)]

### Define df of iHT
iHT.sample.list <- TPMdata[c(1,2,3,4,11,14,15,16,17,18,19,20)]

### Define a dataframe of just adult hypothalamus
aHT.sample.list <- TPMdata[c(7,8,9,10,12)]

### Define a df of iMN
iMN.sample.list <- TPMdata[c(5,6)]

### Define a df of GTEx HT
gtexHT.sample.list <- reference.transcriptomes[122:129]

########################################################################
### Generate lists to select tissues and samples

### Generate a list of filtered reference sets
references.list <- list(ref.full, ref.loose, ref.neutral, ref.tight, ref.supertight)
names(references.list) <- c('Ref.full', 'Ref.loose' ,'Ref.neutral', 'Ref. tight', 'Ref.supertight')

########################################################################
### Plot samples compared against references

### Group of samples against one reference set
spearman.results.list <- list()

### Specify samples
gtexHT.sample.list -> samples
aHT.sample.list -> samples
iHT.sample.list -> samples
iMN.sample.list -> samples

samples <- TPMdata[c(1,2,3)]
samples <- TPMdata[c(4,11,14)]
samples <- TPMdata[c(18,19,20)]


#names(samples.list)
#selection = 10
#samples <- samples.list[[selection]]
#names(samples)[1]
#(target <- targets[selection])
print(str(samples))

### Select a reference set (based on filter level)
names(references.list)
ref.set.num = 1
reference.input <- references.list[[ref.set.num]]              # Specify current reference list and its name
(reference.set.name <- names(references.list)[ref.set.num])    # Declare the name of the rererence filter set

### Define title
(title <- paste(substr(names(samples)[1],1,3), 'vs', reference.set.name))

### Calculate spearman corr. for all samples in group against a chosen reference
spearman.results <- spearman.calc.for.multiple.samples(samples, reference.input)

### Generate a single plot for a single reference set
a <- plot.tissue.match.boxplot(spearman.results, title)

g
a
i
m


### Plot all results on one plot
#spearman.results.list[[length(spearman.results.list)+1]] <- spearman.results

multiboxplot <- plot_grid(g,a,i, labels=c(''), ncol = 3, nrow = 1)
multiboxplot


setwd("z:/Data/Andrew/")
png(filename=paste0('Boxplot',".png"), 
    type="cairo",
    units="in", 
    width=22, 
    height=12, 
    pointsize=12, 
    res=120)
print(multiboxplot)
dev.off()

######################## scratch work

### Reshape data for compatibility with geom_boxplot
new.col <- ncol(spearman.results)+1 
spearman.results[new.col] <- row.names(spearman.results)                    # Add the tissues as a new row
names(spearman.results)[new.col] <- 'ref'                                   # Label the new row
spearman.t <- t(spearman.results[-new.col])                                 # Transpose the data frame, minus new row
spearman.melt <- melt(spearman.t)                                     # Melt transposed data frame
names(spearman.melt) <- c("sample","ref","value")                     # Rename columns of melted data frame
spearman.melt$ref <- factor(spearman.melt$ref, levels = rev(row.names(spearman.results)), ordered = TRUE) # Assign order

### Generate data frame for color and label data
spearman.for.color <- spearman.results                                # Duplicate spearman results
spearman.for.color[-new.col] <- apply(spearman.for.color[-new.col], 1, mean)    # Reassign all values to mean of group
names(spearman.for.color)[2] = 'color.val'
spearman.t <- t(spearman.for.color[-new.col])                               # Transpose again
spearman.melt.color <- melt(spearman.t)                               # Melt again
spearman.melt$color <- spearman.melt.color$value                      # Copy repeating values as 'color column to main df

### Generate data frame for color and label data
spearman.for.color <- spearman.results                                # Duplicate spearman results
spearman.for.color[1:8] <- apply(spearman.for.color[1:8], 1, mean)    # Reassign all values to mean of group
names(spearman.for.color)[2] = 'color.val'
spearman.t <- t(spearman.for.color[-9])                               # Transpose again
spearman.melt.color <- melt(spearman.t)                               # Melt again
spearman.melt$color <- spearman.melt.color$value                      # Copy repeating values as 'color column to main df


g <-   ggplot(data = spearman.melt, aes(x = ref, y = -value, fill = color)) +
  geom_boxplot() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = 'none', panel.background = element_rect(fill = 'grey97'),
        plot.title = element_text(size = 14, face = 'bold', margin = margin(5,0,5,0), hjust = 0.5)) +
  scale_fill_gradientn(colors = c('red','orange','yellow','white')) +
  scale_y_continuous(limits = c(-100,-30), position = 'right') +
  
  labs(title = row.names(spearman.results)[1],
       x = '',
       y = title) +
  geom_text(data = spearman.for.color, aes(x = ref, y = -color.val+8, fill = color.val, label = ref), 
            hjust = 0) +
  coord_flip()
g

g <- ggplot(data = spearman.melt, aes(x = ref, y = -value, fill = color)) +
  geom_boxplot() +
  coord_flip() +
  scale_fill_gradientn(colors = c('red','orange','yellow','white')) +
  scale_y_continuous(limits = c(-100,-30), position = 'right') +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = 'none', panel.background = element_rect(fill = 'grey97'),
        plot.title = element_text(size = 14, face = 'bold', margin = margin(5,0,5,0), hjust = 0.5)) +
  labs(title = row.names(spearman.results)[1],
       x = '',
       y = title) +
  geom_text(data = spearman.for.color, aes(x = ref, y = -color.val+8, fill = color.val, label = ref), 
            hjust = 0)
