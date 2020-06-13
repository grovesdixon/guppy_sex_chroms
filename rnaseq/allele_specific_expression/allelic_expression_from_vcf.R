
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
rm(list=ls())

format_string = "GT:PL:DP:SP:ADF:ADR:AD"


# SELECT DATASET ----------------------------------------------------------

#reticulata
ll = load('rnaseq/vcfs/reticulata_gdat.Rdata');species='reticulata'
sample_cols = colnames(gdat)[grep('SRR', colnames(gdat))]
sexes = read_tsv('metadata/reticulata_rna_runs.txt',
                 col_names = c('run', 'id', 'sex'))
males = sexes %>% 
  filter(sex=='male') %>% 
  pull(run)
females = sexes %>% 
  filter(sex=='female') %>% 
  pull(run)
colnames(gdat)[colnames(gdat) %in% males] = paste(colnames(gdat)[colnames(gdat) %in% males], 'M', sep='_')
colnames(gdat)[colnames(gdat) %in% females] = paste(colnames(gdat)[colnames(gdat) %in% females], 'F', sep='_')
sample_cols = colnames(gdat)[grep('SRR', colnames(gdat))]

#wingei
ll = load('rnaseq/vcfs/wingei_gdat.Rdata')
species='wingei'
sample_cols = colnames(gdat)[grep('Pwin', colnames(gdat))]


#picta
ll = load('rnaseq/vcfs/picta_gdat.Rdata')
species='picta'
sample_cols = colnames(gdat)[grep('P.pic', colnames(gdat))]


# PREP FOR AS EXPRESION ---------------------------------------------------

#now split the genotype data from mpileup to pull out the allele depths
lgdat = gdat %>% 
  dplyr::select(-ID, -QUAL, -FILTER, -INFO, -FORMAT) %>% 
  pivot_longer(cols = all_of(sample_cols),
               names_to = 'indv',
               values_to = 'geno')

sites = paste(lgdat$chr, lgdat$pos, sep='_')
split_geno = lapply(lgdat$geno, function(x) strsplit(x, ':', fixed=TRUE))
ad_string = sapply(split_geno, function(x) return(unlist(x)[7]))
lgdat$ad_string = ad_string
#Darolti et al. compared the 'Major' to 'Minor' allele, so call those as well
dat = lgdat %>% 
  separate(ad_string,
           into = c('d_ref', 'd_alt')) %>% 
  mutate(d_ref = as.numeric(d_ref),
         d_alt = as.numeric(d_alt),
         ref_major = d_ref >= d_alt,
         major = if_else(ref_major,
                         d_ref,
                         d_alt),
         minor = if_else(ref_major,
                         d_alt,
                         d_ref)) %>% 
  dplyr::select(-ref_major)


#assign SDR
sex_chr = 8
left = 0
right = dat %>% 
  filter(chr==sex_chr) %>% 
  pull(pos) %>% 
  max()


#add some additional features
dat = dat %>% 
  mutate(sdr = chr==sex_chr & pos > left & pos < right)
dat$site = paste(dat$chr, dat$pos, sep='_')
dat$sex = if_else(grepl('M',dat$indv),
                  'male',
                  'female')




# TRY TO REPRODUCE FIGURE 3 BASED ON SITES ------------------------------------------

#filter sites where males have reads from both alleles
pmaj_df = dat %>% 
  mutate(sex = if_else(grepl('M', indv),
                       'male',
                       'female'),
         prop_maj = major / (major + minor)) %>% 
  group_by(chr, pos, sdr, gene, sex) %>% 
  summarize(N=n(),
            site_mn_pmaj = mean(prop_maj)) %>% 
  mutate(sex = factor(sex, levels = c('male', 'female')))
  

#combind and plot
pmaj_df %>% 
  ggplot(aes(x=site_mn_pmaj, fill=sdr)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(values = c('grey50', 'goldenrod')) +
  labs(x = 'major allele coverage proportion',
       title = species) +
  facet_wrap(~sex)



# FILTERED VERSION OF FIGURE 3 --------------------------------------------
min_fold_cut = 1

head(dat)
dat$site = paste(dat$chr, dat$pos, sep='_')
dat$sex = if_else(grepl('M',dat$indv),
                  'male',
                  'female')

#filter for sites where all males have at least min_fold_cut
filter_for_min_depth = function(dat, select_sex, min_fold_cut){
  sex_dat = dat %>% 
    filter(sex==select_sex)
  passing_sites = sex_dat %>% 
    group_by(site) %>% 
    summarize(min_count = min(minor)) %>% 
    filter(min_count >= min_fold_cut) %>% 
    pull(site)
  fdat = sex_dat %>% 
    filter(site %in% passing_sites)
}
fm_dat = filter_for_min_depth(dat, 'male', min_fold_cut)
ff_dat = filter_for_min_depth(dat, 'female', min_fold_cut)
fdat = rbind(fm_dat, ff_dat)


#plot by site
fdat %>% 
  mutate(prop_maj = major / (major + minor)) %>% 
  group_by(chr, pos, sdr,sex) %>% 
  summarize(N=n(),
            site_mn_pmaj = mean(prop_maj)) %>% 
  mutate(sex = factor(sex, levels = c('male', 'female'))) %>% 
  ggplot(aes(x=site_mn_pmaj, fill=sdr)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(values = c('grey50', 'goldenrod')) +
  labs(x = 'major allele coverage proportion',
       title = species) +
  facet_wrap(~sex)

#plot for mean gene
fdat %>% 
  mutate(prop_maj = major / (major + minor)) %>% 
  group_by(chr, pos, sdr, gene, sex) %>% 
  summarize(N=n(),
            site_mn_pmaj = mean(prop_maj)) %>% 
  mutate(sex = factor(sex, levels = c('male', 'female'))) %>% 
  ggplot(aes(x=site_mn_pmaj, fill=sdr)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(values = c('grey50', 'goldenrod')) +
  labs(x = 'major allele coverage proportion',
       title = species) +
  facet_wrap(~sex)




# read in xy calls --------------------------------------------------------

read_xy_calls = function(input){
  read_tsv(input) %>% 
    dplyr::select(`#CHROM`, POS, ALT) %>% 
    set_names(c('chr', 'pos', 'hetGamAllele'))
}

add_sex_chrom_calls = function(dat, select_sex, sex_calls){
  sex_dat = dat %>% 
    inner_join(sex_calls, by = c('chr', 'pos')) %>% 
    mutate(hetGamDepth = if_else(hetGamAllele==ref,
                                 d_ref,
                                 d_alt),
           homoGamDepth = if_else(hetGamAllele==ref,
                                  d_alt,
                                  d_ref),
           propHomoGam = homoGamDepth / (homoGamDepth+hetGamDepth)) %>% 
    filter(sex==select_sex)
  return(sex_dat)
}

xy_calls_in = paste('rnaseq/callXY/', species, '_hcut8_y.vcf', sep='')
zw_calls_in = paste('rnaseq/callXY/', species, '_hcut8_w.vcf', sep='')
xy_calls = read_xy_calls(xy_calls_in)
zw_calls = read_xy_calls(zw_calls_in)
xy_dat = add_sex_chrom_calls(dat, 'male', xy_calls)
zw_dat = add_sex_chrom_calls(dat, 'female', zw_calls)

#plot density
rbind(xy_dat, zw_dat) %>% 
  group_by(sdr, sex, gene) %>% 
  summarize(N=n(),
            site_mn_ratio = median(propHomoGam)) %>% 
  ggplot(aes(x=site_mn_ratio, fill=sdr)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(values = c('grey50', 'goldenrod')) +
  labs(x = 'major allele coverage proportion',
       title = species) +
  lims(x=c(0,1)) +
  geom_vline(xintercept = 0.5, lty=2) +
  facet_wrap(~sex)


#plot along sex chromosome
xy_dat %>% 
  group_by(chr, sdr, sex, gene) %>% 
  summarize(N=n(),
            site_mn_ratio = median(propHomoGam),
            med_pos = median(pos)) %>% 
  filter(chr==8) %>% 
  ggplot(aes(x=med_pos, y=site_mn_ratio)) +
  geom_point() +
  geom_smooth()


