#!/usr/bin/env Rscript
#call_xy.R

library(tidyverse)
library(optparse)
options(stringsAsFactors = FALSE) #important for WGCNA

option_list = list(
  
  make_option(c("--vcf"), type="character", 
              help="input vcf file from a single chromosome with all depth (DP) in same place in FORMAT column"),
  make_option(c("--window_size"), type="numeric",  default = 10000,
              help="size of the windows in bp. Default is 10Kb"),
  make_option(c("--o"), type="character", 
              help="output file name")
  )

print("Parsing arugments...")
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
vcf_input = opt$vcf
window_size = opt$window_size
outfile_name = opt$o


# load data ---------------------------------------------------------------

print(paste('Reading in ', vcf_input, '...', sep=''))
vcf_dat = read_tsv(vcf_input, comment = '#')
geno_dat = vcf_dat %>% 
  select(10:ncol(vcf_dat))


# find the AD in FORMAT ---------------------------------------------------

get_ad_index = function(x){
  vec = unlist(strsplit(x, ':'))
  grep('^DP$', vec)
}

depth_indices = sapply(vcf_dat$FORMAT, function(x) get_ad_index(x))
di = unique(depth_indices)
if (length(di) > 1){
  print('ERROR. AD must be in same spot for all variants for this script to work.')
  exit()
}

#also get the chromosome
vcf_chr = unique(vcf_dat$CHROM)
if (length(vcf_chr) > 1){
  print('ERROR. The VCF is expected to be from a single chromosome.')
  exit()
}

# sum the allele depths ---------------------------------------------------
print('Pulling out depths...')
depth_dat = vcf_dat %>% 
  select(1:9) %>% 
  data.frame()
for (cn in colnames(geno_dat)){
  v=geno_dat %>% 
    pull(cn)
  depths = sapply(v, function(x) return(strsplit(x, ':')[[1]][di])) %>% 
    as.numeric()
  depth_dat[,cn] = depths
}


# get window means --------------------------------------------------------

print('Getting depths for means...')
#assign windows
ml = max(vcf_dat$POS)
bounds = seq(0, ml, by = window_size)
long_depth = depth_dat %>% 
  pivot_longer(10:ncol(depth_dat),
               names_to = 'sample',
               values_to = 'depth')
bins = cut(long_depth$POS, breaks = bounds)
long_depth$window = as.character(bins)


#get median depth by window
library(plotrix)
wdat = long_depth %>% 
  group_by(sample, window) %>% 
  summarize(window_mean = mean(depth, na.rm=TRUE),
            window_median = median(depth, na.rm=TRUE),
            window_sd = sd(depth, na.rm = TRUE),
            window_sterr = std.error(depth),
            nsites = n()) %>% 
  mutate(chr = vcf_chr) %>% 
  as_tibble



# write out ---------------------------------------------------------------
print(paste('Writing out results to ', outfile_name, '...', sep=''))
wdat %>% 
  write_tsv(path = outfile_name)