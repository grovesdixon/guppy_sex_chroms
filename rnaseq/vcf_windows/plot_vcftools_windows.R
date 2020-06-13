#plot_vcftools_windows.R
library(tidyverse)
library(plotrix)
library(cowplot)
theme_set(theme_cowplot())
rm(list=ls())



# UPLOAD VCF WINDOW WRAPPER RESULTS ---------------------------------------

#loop through the species
species = c('reticulata',
            'wingei',
            'picta',
            'latipinna',
            'Gambusia')

sp_datlist = list()
for (sp in species){
  print(sp)
  infile = paste('rnaseq/vcf_windows/',sp, '_100Kb_win_res.tsv', sep='')
  dat = read_tsv(infile) %>% 
    dplyr::select(-X1) %>% 
    mutate(species=sp,
           mb = BIN_START/1e6,
           mf_pi_ratio = log( (PI_male / PI_female), 2),
           mf_snpCount_ratio = log( (SNP_COUNT_male / SNP_COUNT_female), 2)) %>% 
    select(CHROM,
           mb,
           N_VARIANTS,
           pMaleSpecific,
           pFemaleSpecific,
           MEAN_FST,
           meanDiffA1,
           mf_pi_ratio,
           mf_snpCount_ratio,
           species)
  sp_datlist[[sp]] = dat
}

dat = purrr::reduce(sp_datlist, rbind)
dat


get_95_cut = function(df, col){
  df=data.frame(df)
  v=df[,col]
  q = quantile(v, probs=seq(0,1,0.05), na.rm=TRUE)
  cut = q['95%']
  return(v>cut)
}

cols = c('N_VARIANTS',
         'pMaleSpecific',
         'pFemaleSpecific',
         'MEAN_FST',
         'meanDiffA1',
         'mf_pi_ratio',
         'mf_snpCount_ratio')

add_cuts = function(df){
  mod_df = df
  for (c in cols){
    print(c)
    mod_df[,paste(c,'95',sep='.')] = get_95_cut(df, c)
  }
  return(mod_df)
}

wcut_list = lapply(sp_datlist, function(x) add_cuts(x))

dat = purrr::reduce(wcut_list, rbind)
dat

dat$species = factor(dat$species, levels = species)

# PLOT STATS --------------------------------------------------------------

plot_sex_diff_col = function(df, col_to_plot, CHR=8, YLIM=FALSE){
  v=df[,col_to_plot]
  q = quantile(v, probs=seq(0,1,0.05), na.rm=TRUE)
  cut = q['95%']
  cut_col = paste(col_to_plot, '95', sep='.')
  plt = df %>% 
    filter(CHROM==CHR) %>% 
    mutate(revMb = max(mb)-mb) %>% 
    ggplot(aes_string(x='revMb', y=col_to_plot)) +
    geom_point(aes_string(color=cut_col)) +
    geom_smooth(se=FALSE) +
    scale_color_manual(values=c('black', 'red')) +
    facet_wrap(~species, nrow=length(species)) 
  if (length(YLIM)>1){
    plt = plt + lims(y=YLIM)
  }
  return(plt)
}

CHR=8
plot_sex_diff_col(dat, 'N_VARIANTS', CHR=CHR)
plot_sex_diff_col(dat, 'pMaleSpecific', CHR=CHR, YLIM=c(0,0.08))
plot_sex_diff_col(dat, 'pFemaleSpecific', CHR=CHR, YLIM=c(0,0.08))
plot_sex_diff_col(dat, 'MEAN_FST', CHR=CHR)
plot_sex_diff_col(dat, 'meanDiffA1', CHR=CHR)
plot_sex_diff_col(dat, 'mf_pi_ratio', CHR=CHR)
plot_sex_diff_col(dat, 'mf_snpCount_ratio', CHR=CHR)

vcf_dat = dat
save(vcf_dat, file='figure_plotting/vcf_dat.Rdata')



# PLOT TO MATCH DAROLTI ---------------------------------------------------

spp = 'picta'


#plot all along the sex chromosome
dat %>% 
  filter(species==spp) %>% 
  select(c('CHROM', 'mb', 'pMaleSpecific', 'pFemaleSpecific', 'MEAN_FST', 'mf_pi_ratio')) %>% 
  pivot_longer(pMaleSpecific:mf_pi_ratio,
               names_to = 'stat',
               values_to = 'value') %>% 
  filter(CHROM==8) %>% 
  mutate(stat = factor(stat, levels=c('pMaleSpecific', 'pFemaleSpecific', 'MEAN_FST', 'mf_pi_ratio'))) %>% 
  ggplot(aes(x=mb, y=value)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  facet_wrap(~stat, scales = 'free_y', nrow=4)


#barplot
chrdat %>% 
  pivot_longer(malePi:pMaleSpecific,
               names_to = 'stat',
               values_to = 'value') %>% 
  mutate(stat = factor(stat, levels=c('malePi', 'femalePi', 'mf_pi_ratio', 'fst', 'pMaleSpecific'))) %>% 
  ggplot(aes(x=CHROM, y=value, fill=CHROM)) +
  geom_bar(stat='identity') +
  facet_wrap(~stat, scales='free')


#plot along chromosomes
ycol = 'mf_pi_ratio'
dat %>% 
  filter(species==spp) %>% 
  ggplot(aes_string(x='mb', y=ycol)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  geom_hline(yintercept = 0, lty=2, color='red') +
  labs(y=bquote(F[ST]~'M:F'), x='position (Mb)') +
  facet_wrap(~CHROM, scales = 'free_x')



#density plot of male:female pi ratio
sex_chrom = dat %>% 
  filter(species==spp,
         CHROM==8) %>% 
  ggplot(aes_string(x='mf_pi_ratio')) +
  geom_density() +
  geom_vline(xintercept = c(-1, 0))



#density plot of male:female pi ratio
sex_chrom = dat %>% 
  filter(species==spp) %>% 
  ggplot(aes_string(x='mf_pi_ratio')) +
  geom_density() +
  geom_vline(xintercept = 0, lty=2) +
  geom_vline(xintercept = -1, lty=1) +
  labs(x=bquote("log"[2]~'Male:Female Pi Ratio')) +
  facet_wrap(~CHROM)



#look at heterozygosity
het_infile = paste('rnaseq/callXY/', spp, '_heterozygosities.tsv', sep='')
hdat = read_tsv(het_infile) %>% 
  mutate(sex = if_else(grepl('M', INDV),
                       'male',
                       'female'),
         ho = 1-(`O(HOM)` / N_SITES))

#overall
hdat %>% 
  ggplot(aes(x=sex, y=ho)) +
  geom_boxplot()

#inbreeding coefficient
hdat %>% 
  ggplot(aes(x=sex, y=F)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty=2)


#by individual and chromosome
inds = hdat %>% 
  arrange(sex) %>% 
  pull(INDV) %>% 
  unique()

hdat %>% 
  mutate(ind = factor(INDV, levels = inds),
         chr = factor(chr, levels = paste('chr', 1:24, sep=''))) %>% 
  ggplot(aes(x=ind, fill=sex, y = ho)) +
  geom_bar(stat='identity') +
  facet_wrap(~chr)
