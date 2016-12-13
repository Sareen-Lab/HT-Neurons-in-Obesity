######################################################################################
###### Proteomic analysis for networks, 2016-02-17
######################################################################################
### Functions

convert.char.column.to.num <- function(dataframe,column) {
  new.vector <- c()
  for (value in dataframe[,column]) {
    new.value <- as.numeric(strsplit(as.character(value),"/")[[1]][1])
    new.vector <- c(new.vector,new.value)
  }
  dataframe[column] <- new.vector
  return(dataframe)
}
######################################################################################
### Upload data

setwd(dir = "Z:/Uthra/HT paper/Bioinformatics figures/IPA analysis/network_files_for_cytoscape/")
networks <- read.csv("z:/Uthra/HT paper/Bioinformatics figures/IPA analysis/network_files_for_cytoscape/Pathway_table-iHT--Obs-v-Ctr--8000_nc_ETROG.csv")

networks <- networks[-4]
length <- nrow(networks)

######################################################################################
### Format new columns

for(column in 4:7){
  networks <- convert.char.column.to.num(networks,column)
}

######################################################################################
### Calculate shared genes

outputNetwork <- data.frame(sources=character(),target=character(),
                     interaction=character(),boolean=character(),
                     string=character(),shared=double(),
                     stringsAsFactors=FALSE)

for (i in 1:(length-1)) {                                           # Run through each node in turn
  sourceName <- toString(networks$Ingenuity.Canonical.Pathways[i])        # Define current source
  molecules <- networks$Molecules[i]                            # Generate a string of molecules in the current node
  moleculesA <- unlist(strsplit(toString(molecules),","))       # Convert string to list
  for (j in 1:(length-i)) {                                     # Run through each of the nodes following the current one
    targetName <- toString(networks$Ingenuity.Canonical.Pathways[i+j])      # Define current source
    molecules <- networks$Molecules[i+j]                        # Generate a string of molecules
    moleculesB <- unlist(strsplit(toString(molecules),","))     # Convert string to list
    count <- length(intersect(moleculesA,moleculesB))           # Calculate the number of genes shared between the nodes
    row <- c(sourceName,targetName,"cooccurrence","TRUE","ABC",count)
    outputNetwork[length(outputNetwork[,1])+1,] <- row
  }
}

outputNode <- networks[1:7]
names(outputNode) <- c("sources","pVal","ratio","Downregulated","No.change","Upregulated","No.overlap")

######################################################################################
### Output shared genes

setwd(dir = "Z:/Uthra/HT paper/Bioinformatics figures/IPA analysis/network_files_for_cytoscape/")
names(outputNetwork) <- c("source","target","interaction","boolean attribute","string attribute","floating point attribute")
write.table(outputNetwork,paste0("NET--iHT-8000-for_cytoscape_",substr(weekdays(Sys.Date()),1,4),"-",format(Sys.Date(),"%b-%d"),".txt"),row.names=FALSE,sep="\t",quote=FALSE)

write.table(outputNode,paste0("NODE--iHT-8000-for_cytoscape_",substr(weekdays(Sys.Date()),1,4),"-",format(Sys.Date(),"%b-%d"),".txt"),row.names=FALSE,sep="\t",quote=FALSE)


