#!/usr/bin/env Rscript
#pca_chrom_snps.R
#Groves Dixon
#9-11-17
#last updated 4-22-20 -- added step to handle of missing data
#Use Adegenet to build PCAs from a set of vcf files

#requires file with list of file prefixes
#library(plotrix)
library(vcfR)
library(adegenet)


#parse arguments
args = commandArgs(trailingOnly=TRUE)
vcfInput = args[1] #table with sample labels
N.CORES = as.numeric(args[2])


print("Running R script pca_chrom_snps.R...")
print(paste("Number of cores to use =", N.CORES))


#upload the vcf
print(paste("Loading VCF file", vcfInput))
gll=vcfR2genlight(read.vcfR(vcfInput))
print("Done loading VCF.")
print(gll)

#check for missing data
missing = glNA(gll)
has_missing = missing > 0
n_has_missing = sum(has_missing)
tot_loci = ncol(gll)
pct_missing = round((n_has_missing / tot_loci), 3)*100
pct_string = paste('(', pct_missing, '%)', sep='')
good = which(!has_missing)
if (n_has_missing > 0){
	print(paste('WARNING.', n_has_missing, 'of', tot_loci, pct_string, 'loci have at least one NA and will be removed before PCA'))
	gll = gll[,good]
	print('Summary of new genlight object:')
	print(gll)
}


#run PCA
print("Running PCA...")
pca=glPca(gll,nf=6, n.cores=N.CORES)
p=pca$scores
print("Done.")

#output
outFile = paste(vcfInput, "pca.Rdata", sep="_")
print(paste("Saving results as", outFile))
save(pca, file=outFile)