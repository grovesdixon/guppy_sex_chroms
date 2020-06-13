#Poecilia reticulata

#set variables
TITLE="P.picta : G. holbrooki"
SPECIES=c("Gambusia_holbrooki",
          "Micropoecilia_picta")
PREFIX="Ppicta_v_Gholbrooki"
CONTRAST = c('Organism', 'Micropoecilia_picta', 'Gambusia_holbrooki')


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

#make sex subsets
fcoldata = coldata[coldata$sex=='female',]
mcoldata = coldata[coldata$sex=='male',]
fcounts = counts[,rownames(fcoldata)]
mcounts = counts[,rownames(mcoldata)]
