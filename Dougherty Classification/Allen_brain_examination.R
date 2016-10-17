### Allen Brain Atlas file check -- Andrew R Gross -- 2016/10/11
### Examining files for hypothalamus


########################################################################
### Header
########################################################################


########################################################################
### Functions
########################################################################



########################################################################
### Input
########################################################################

setwd('c://Users/grossar/Bioinform/DATA/Allen Brain Atlas/')

### Human Brain - RNAseq
human.brain.rna.61 <- read.csv('Human Brain/RNAseq datasets/rnaseq_donor9861/RNAseqTPM.csv', row.names = 1, header = FALSE)
samples.61 <- read.csv('Human Brain/RNAseq datasets/rnaseq_donor9861/SampleAnnot.csv')
ontology <- read.csv('Human Brain/RNAseq datasets/rnaseq_donor9861/Ontology.csv')

present.ontologies <- ontology[match(samples.61$ontology_structure_id,ontology$id),]

human.brain.rna.21 <- read.csv('Human Brain/RNAseq datasets/rnaseq_donor10021/RNAseqTPM.csv', row.names = 1, header = FALSE)
samples.21 <- read.csv('Human Brain/RNAseq datasets/rnaseq_donor10021/SampleAnnot.csv')

present.ontologies <- ontology[match(samples.21$ontology_structure_id,ontology$id),]

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


### Developing Brain

########################################################################
### Format
########################################################################
