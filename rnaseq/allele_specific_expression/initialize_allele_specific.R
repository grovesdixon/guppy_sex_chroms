#initialize_counts.R
#prepare reads for DESeq and get gist of data

library('DESeq2')
library('cowplot')
theme_set(theme_cowplot())
library('tidyverse')
rm(list=ls())
source('rnaseq/scripts/rnaseq_functions.R')


#upload read counts from htseq
counts = read.table("rnaseq/allele_specific_expression/feature_counts_out.tsv", header = T, row.names='Geneid')
counts = counts[,6:ncol(counts)]
head(counts)
length(unique(colnames(counts)))


# SEPARATE OUT THE DATA BY SPECIES ----------------------------------------

cols = colnames(counts)

############ picta
picta = counts %>% 
  dplyr::select(contains('picta'))
cols = colnames(picta)
mod_cols = sub('.bam', '', cols, fixed=TRUE)
mod_cols = sub('T_sortedByName.', '', mod_cols, fixed=TRUE)
mod_cols = sub('_genome1', '_X', mod_cols, fixed=TRUE)
mod_cols = sub('_genome2', '_Y', mod_cols, fixed=TRUE)
mod_cols = sub('picta_', '', mod_cols, fixed=TRUE)
samples = sapply(mod_cols, function(x) strsplit(x, '_')[[1]][1])
table(samples)

spp_counts = picta
colnames(picta) = mod_cols

sum_paired_and_singletons = function(spp_counts, samples){
  summed_counts = data.frame(row.names = rownames(spp_counts))
  for (s in unique(samples)){
    print(s)
    paired_x = paste(s, 'X', sep='_')
    singleton_x = paste(s, 'X_st', sep='_')
    paired_y = paste(s, 'Y', sep='_')
    singleton_y = paste(s, 'Y_st', sep='_')
    summed_counts[,paste(s, 'X', sep='_')] = spp_counts[,paired_x] + spp_counts[,singleton_x]
    summed_counts[,paste(s, 'Y', sep='_')] = spp_counts[,paired_y] + spp_counts[,singleton_y]
  }
  return(summed_counts)
}

build_coldata = function(dat){
  samples = colnames(dat)
  individual = sapply(samples, function(x) strsplit(x, '_')[[1]][1])
  sex = if_else(grepl('M', samples),
                'male',
                'female')
  chrom = if_else(grepl('_Y', samples),
                  'Y',
                  'X')
  coldata = data.frame(individual,
                       sex,
                       chrom,
                       row.names = samples)
  return(coldata)
}

summed_counts = sum_paired_and_singletons(picta, samples)
coldata = build_coldata(summed_counts)
sum(rownames(coldata) == colnames(summed_counts))==ncol(summed_counts)
save(summed_counts, coldata, file='rnaseq/allele_specific_expression/picta_xy_counts.Rdata')

############ wingei
wingei = counts %>% 
  dplyr::select(contains('wingei'))
cols = colnames(wingei)
mod_cols = sub('.bam', '', cols, fixed=TRUE)
mod_cols = sub('sortedByName.', '', mod_cols, fixed=TRUE)
mod_cols = sub('_genome1', '_X', mod_cols, fixed=TRUE)
mod_cols = sub('_genome2', '_Y', mod_cols, fixed=TRUE)
mod_cols = sub('wingei_', '', mod_cols, fixed=TRUE)
samples = sapply(mod_cols, function(x) strsplit(x, '_')[[1]][1])
table(samples)

colnames(wingei) = mod_cols
summed_counts = sum_paired_and_singletons(wingei, samples)

coldata = build_coldata(summed_counts)
sum(rownames(coldata) == colnames(summed_counts))==ncol(summed_counts)
save(summed_counts, coldata, file='rnaseq/allele_specific_expression/wingei_xy_counts.Rdata')

############ reticulata
reticulata = counts %>% 
  dplyr::select(contains('reticulata'))
cols = colnames(reticulata)
mod_cols = sub('.bam', '', cols, fixed=TRUE)
mod_cols = sub('lab_RNA.sortedByName.', '', mod_cols, fixed=TRUE)
mod_cols = sub('_genome1', '_X', mod_cols, fixed=TRUE)
mod_cols = sub('_genome2', '_Y', mod_cols, fixed=TRUE)
mod_cols = sub('reticulata_', '', mod_cols, fixed=TRUE)
mod_cols = sub('F_', 'F', mod_cols, fixed=TRUE)
mod_cols = sub('M_', 'M', mod_cols, fixed=TRUE)
samples = sapply(mod_cols, function(x) strsplit(x, '_')[[1]][1])
table(samples)
colnames(reticulata) = mod_cols
summed_counts = sum_paired_and_singletons(reticulata, samples)

coldata = build_coldata(summed_counts)
sum(rownames(coldata) == colnames(summed_counts))==ncol(summed_counts)
save(summed_counts, coldata, file='rnaseq/allele_specific_expression/reticulata_xy_counts.Rdata')




