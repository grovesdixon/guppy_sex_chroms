library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
rm(list=ls())

#choose teh input vcf files
input_files = c('rnaseq/vcfs/reticulata_genome_maf.recode.vcf',
                'rnaseq/vcfs/wingei_genome_maf.recode.vcf',
                'rnaseq/vcfs/picta_genome_maf.recode.vcf')
dat_names = c('reticulata',
              'wingei',
              'picta')
names(input_files) = dat_names


#read the vcfs into list
load_vcf = function(input_vcf){
  vdat = read_tsv(input_vcf,
                  comment = '##') %>% 
    rename(chr = `#CHROM`,
           pos = POS,
           ref = REF,
           alt = ALT)
  return(vdat)
}

vcf_datlist = map(input_files, load_vcf)


#format and add gene data and save it
for (n in dat_names){
  print('=========')
  print(n)
  vdat = vcf_datlist[[n]]
  #grab sample columns and format string
  sample_cols = colnames(vdat)[10:ncol(vdat)]
  format_string = unique(vdat$FORMAT)
  if (length(format_string) > 1){
    print('ERROR! There is more than one type of FORMAT string in this VCF.')
    print('This script only works for VCFs with a single FORMAT string for all variants')
  }
  
  #get the index for the allele depth (AD in the FORMAT string)
  #note that this script only works when all rows in the VCF have the same FORMAT string
  f = vdat$FORMAT
  splits = lapply(f, function(x) strsplit(x, ':', fixed=TRUE))
  get_index = function(x){
    s = unlist(x)
    is = 1:length(s)
    i = is[s=='AD']
    return(i)
  }
  vdat$ad_index = sapply(splits, function(x) get_index(x))
  ad_indices = sapply(splits, function(x) get_index(x))
  ad_index = unique(ad_indices)
  
  
  #assign the genes
  gene_coords = read_tsv('metadata/geneCoords.tsv',
                         col_names = c('chr', 'start', 'end', 'gene')) 
  gdat = data.frame()
  for (i in 1:nrow(gene_coords)){
    if (i %% 1000 == 0){
      print(paste(i, 'of', nrow(gene_coords), 'complete'))
    }
    row = gene_coords[i,]
    chr = row$chr
    vsub = vdat %>% 
      filter(chr == row$chr,
             pos >= row$start,
             pos <= row$end) %>% 
      mutate(gene = row$gene) %>% 
      data.frame()
    gdat = rbind(gdat, vsub)
  }
  outname = paste('rnaseq/vcfs/', n, '_gdat.Rdata', sep='')
  save(gdat, sample_cols, format_string, file=outname)
}
