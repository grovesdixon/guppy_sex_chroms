#Poecilia reticulata

#set variables
TITLE='Poecilia reticulata'
SPECIES=c('Poecilia_reticulata')
PREFIX='PreticulataOnly'


#load all RNA data
ll = load('rnaseq/deseqInput.Rdata')
ll
coldata
head(counts)
sum(colnames(counts) == rownames(coldata))==ncol(counts)


#make subsets
keep = rownames(coldata)[coldata$Organism %in% SPECIES]
length(keep)
coldata = coldata[keep,]
counts=counts[,keep]
dim(coldata)
dim(counts)
sum(colnames(counts) == rownames(coldata))==ncol(counts)