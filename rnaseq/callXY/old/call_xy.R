#call_xy.R

library(tidyverse)
library(optparse)
options(stringsAsFactors = FALSE) #important for WGCNA

option_list = list(
  
  make_option(c("--m_prefix"), type="character", 
              help="prefix for male *afs.tsv and *hwe.tsv files"),
  make_option(c("--f_prefix"), type="character", 
              help="prefix for female *afs.tsv and *hwe.tsv files"),
  make_option(c("--out_prefix"), type="character",
              help="prefix output files"),
  make_option(c("--hcut"), type="numeric", default = 0.51,
              help="heterozygosity cutoff. The heterogametic sex must have heterozygosity >= to this value for the allele to be considered W or Y-linked"),
  make_option(c("--hetgam_freq_cut"), type="numeric", default = 0.3,
              help="frequency cutoff for heterogametic sex. Frequency must be >= to this value in the heterogametic sex for the allele to be considered W or Y-linked"),
  make_option(c("--homgam_freq_cut"), type="numeric", default = 0.0,
              help="frequency cutoff for homogametic sex. Frequency must be <= to this value in the homogametic sex for the allele to be considered W or Y-linked"),
  
)

print("Parsing arugments...")
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
m_prefix = opt$m_prefix
f_prefix = opt$f_prefix
out_prefix = opt$out_prefix
het_cut = opt$hcut
heterogametic_freq_cut = opt$hetgam_freq_cut
homogametic_freq_cut = opt$homgam_freq_cut


read_freq = function(infile){
  d = read_tsv(infile, col_names = c('chr', 'pos', 'Nallele', 'Nchr', 'a1', 'a2'))
  ref_split = sapply(d$a1, function(x) strsplit(x, ':', fixed=TRUE))
  alt_split = sapply(d$a2, function(x) strsplit(x, ':', fixed = TRUE))
  ref_allele = sapply(ref_split, function(x) x[1])
  alt_allele = sapply(alt_split, function(x) x[1])
  ref_frq = sapply(ref_split, function(x) x[2]) %>% 
    as.numeric()
  alt_frq = sapply(alt_split, function(x) x[2]) %>% 
    as.numeric()
  d %>% 
    mutate(ref = ref_allele,
           alt = alt_allele,
           rfrq = ref_frq,
           afrq = alt_frq,
           site = paste(chr, pos, sep='-')) %>% 
    select(chr, pos, ref, alt, rfrq, afrq, site)
}

mfrq = read_freq(paste(m_prefix, '_afs.tsv', sep=''))
ffrq = read_freq(paste(f_prefix, '_afs.tsv', sep=''))
print(dim(mfrq))
print(dim(ffrq))

#function to read in the hwe data
read_hwe = function(infile){
  dat = read_tsv(infile, col_names = c('chr', 'pos', 'obs')) %>% 
    separate(obs,
             into = c('Nhomo1', 'Nhet', 'Nhomo2'),
             convert = TRUE) %>% 
    mutate(Nind = Nhomo1+Nhet+Nhomo2,
           het = Nhet/Nind) %>% 
    select(chr, pos, Nhet, het)
}

mhw = read_hwe(paste(m_prefix, '_hwe.tsv', sep=''))
fhw = read_hwe(paste(f_prefix, '_hwe.tsv', sep=''))
print(dim(mhw))
print(dim(fhw))



print('checking sites line up...')
check1 = sum(mfrq$site==ffrq$site)==nrow(mfrq)
check2 = sum(mfrq$ref==ffrq$ref)==nrow(ffrq)
check3 = sum(mhw$pos == fhw$pos) == nrow(fhw)
check4 = sum(mhw$pos == mfrq$pos) == nrow(mfrq)
check = check1 & check2 & check3 & check4
if (check){
  print('sites match up. Continuing to call XY...')
} else{
  print("Error! Sites don't match.")
  print('Exiting.')
  quit()
}

#combine into single df
dat = mfrq %>% 
  select(-site) %>% 
  rename(m_ref_frq = rfrq,
         m_alt_frq = afrq) %>% 
  data.frame()
dat$f_ref_frq = ffrq$rfrq
dat$f_alt_frq = ffrq$afrq
dat$m_het = mhw$het
dat$f_het = fhw$het
head(dat)
dim(dat)

hetgam_prefix = 'm'
homogam_prefix = 'f'

get_ylike = function(dat, hetgam_prefix, homogam_prefix, y_freq_cut, x_freq_cut){
  refY = paste(hetgam_prefix, 'ref_frq', sep='_')
  refX = paste(homogam_prefix, 'ref_frq', sep='_')
  altY = paste(hetgam_prefix, 'alt_frq', sep='_')
  altX = paste(homogam_prefix, 'alt_frq', sep='_')
  hetY = paste(hetgam_prefix, 'het', sep='_')
  hetX = paste(homogam_prefix, 'het', sep='_')
  pass_hetY = dat[,hetY] >= het_cut
  pass_hetX = dat[,hetX] >= het_cut
  pass_refY = dat[,refY] >= y_freq_cut
  pass_refX = dat[,refX] <= x_freq_cut
  pass_ref = pass_refY & pass_refX
  pass_altY = dat[,altY] >= y_freq_cut
  pass_altX = dat[,altX] <= x_freq_cut
  pass_alt = pass_altY & pass_altX
  full_pass_ref = pass_ref & pass_hetY
  full_pass_alt = pass_alt & pass_hetY
  full_pass_ref[is.na(full_pass_ref)] <- FALSE
  full_pass_alt[is.na(full_pass_alt)] <- FALSE
  ylike = full_pass_ref | full_pass_alt
  dat$Yallele = if_else(full_pass_ref,
                        'ref',
                        'alt')
  dat$YalleleCheck = if_else(full_pass_alt,
                             'alt',
                             'ref')
  
  
  
  print(paste('total SNPs =', nrow(dat)))
  print(paste('total SNPs passing heterozygosity rate filter for heterogametic sex', sum(pass_hetY, na.rm=TRUE)))
  print(paste('total SNPs passing heterozygosity rate filter for homogametic sex', sum(pass_hetX, na.rm=TRUE)))
  print(paste('total SNPs passing reference allele frequency rate cuts', sum(pass_ref, na.rm=TRUE)))
  print(paste('total SNPs passing alternative allele frequency rate cuts', sum(pass_alt, na.rm=TRUE)))
  print(paste('total reference alleles passing frequency AND het cuts', sum(full_pass_ref, na.rm=TRUE)))
  print(paste('total alternative alleles passing frequency AND het cuts', sum(full_pass_alt, na.rm=TRUE)))
  print(paste('total Y-like SNPs =', sum(ylike, na.rm=TRUE)))
  ydat = dat[ylike,] %>% 
    as_tibble()
  return(ydat)
}


ydat = get_ylike(dat=dat,
                 hetgam_prefix = 'm',
                 homogam_prefix = 'f',
                 y_freq_cut = heterogametic_freq_cut,
                 x_freq_cut = homogametic_freq_cut)

wdat = get_ylike(dat=dat,
                 hetgam_prefix = 'f',
                 homogam_prefix = 'm',
                 y_freq_cut = heterogametic_freq_cut,
                 x_freq_cut = homogametic_freq_cut)

dim(ydat)
dim(wdat)


# OUTPUT SUMMARY ----------------------------------------------------------

summary_outname = paste(out_prefix,'ylike_summary.tsv', sep='_')

ydat$source = m_prefix
wdat$source = f_prefix
ylike_dat = rbind(ydat, wdat)
ylike_sum = ylike_dat %>% 
  group_by(source, chr) %>% 
  summarize(N_ylike = n())
ylike_sum %>% 
  write_tsv(path=summary_outname)




# OUTPUT SNPS FOR SNP SPLIT -----------------------------------------------

#function to output a snp file in format for running SNP split 
#(https://github.com/FelixKrueger/SNPsplit/blob/master/SNPsplit_User_Guide.md)
#EG:
# ID	   Chr   Position  SNP value	Ref/SNP
# 18819008	 5	48794752	  1	      C/T 
# 40491905	11	63643453	  1	      A/G 
# 44326884	12	96627819	  1	      T/A 

#note here that the Ref indicates the allele for one genome,
#and SNP representst the allele for the other genome for allele-specific expression.
#Hence for an X and a Y, put the X allele as Ref and the Y-allele as SNP
#regardless of the reference/alternative state from the VCF

make_snp_file = function(sex_dat){
  ydat %>% 
    mutate(ID=paste(chr, pos, sep='_'),
           Chr = chr,
           Position = pos,
           `SNP value` = 1,
           `Ref/SNP` = if_else(Yallele=='alt',
                               paste(ref,alt,sep='/'),
                               paste(alt,ref,sep='/'))) %>% 
    dplyr::select(ID, Chr, Position, `SNP value`, `Ref/SNP`)
}

y_snps = make_snp_file(ydat)
w_snps = make_snp_file(wdat)

y_snps %>% 
  write_tsv(path=paste(out_prefix, 'y_snpfile.tsv', sep='_'))

w_snps %>% 
  write_tsv(path=paste(out_prefix, 'w_snpfile.tsv', sep='_'))





