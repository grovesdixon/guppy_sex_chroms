#plot_figure1.R
rm(list=ls())
source('functions.R')
species_names

ALPHA = 0.75

# LOAD FOLD COVERAGE DATA -------------------------------------------------
#build this with dna/fold_coverage.R

ll=load('dna/fold_coverage/full_maculatus/fold_dat.Rdata') #from full mapping to maculatus reference
#other options that were checked, but not used in final figure:
#ll=load('dna/fold_coverage/full_hellerii/fold_dat.Rdata') #from full mapping to maculatus reference
#ll=load('dna/fold_coverage/quick_maculatus/fold_dat.Rdata') #from 1 million reads mapped to hellerii reference
#ll=load('dna/fold_coverage/quick_hellerii/fold_dat.Rdata') #from 1 million reads mapped to hellerii reference

ll
fold_dat = fold_dat %>% 
  mutate(species = factor(species, levels=species_names),
         mb=start/1e6)


# LOAD W-LIKE FREQUENCY DATA ----------------------------------------------
#load results based on DNA (run with the 'both' pipeline)
ll=load('both_dna_rna/callXY/wdens_dat.Rdata')
ll
wdens_dat = ydens_dat %>% 
  mutate(species = factor(species, levels = species_names),
         mb = starts/1e6,
         wlength = ends-starts,
         ylike_density = nYlike/wlength)
remove(ydens_dat)


# LOAD Y-LIKE FREQUENCY DATA ----------------------------------------------
#build these with y_like_density.R

# #load results based on RNA (kept here for reference)
# ll=load('rnaseq/callXY/ydens_dat.Rdata')
# ll
# ydens_dat = ydens_dat %>% 
#   mutate(species = factor(species, levels = species_names),
#          mb = starts/1e6,
#          wlength = ends-starts,
#          ylike_density = nYlike/wlength)


#load results based on DNA (run with the 'both' pipeline)
ll=load('both_dna_rna/callXY/ydens_dat.Rdata')
ll
ydens_dat = ydens_dat %>% 
  mutate(species = factor(species, levels = species_names),
         mb = starts/1e6,
         wlength = ends-starts,
         ylike_density = nYlike/wlength)


# LOAD VCF DATA -----------------------------------------------------------

ll=load('figure_plotting/dna_vcf_dat.Rdata') #build this with rnaseq/plot_vcftools_windows.R

vcf_dat = vcf_dat %>% 
  mutate(species = factor(species, levels=species_names))


# LOAD BASIC GE DATA ------------------------------------------------------

ll=load('figure_plotting/ge_res.Rdata') #build this with rnaseq/chromosome_comparisons.R
ll
ge_res = ge_res %>% 
  mutate(mb=start/1e6,
         species = factor(species, levels=species_names))


# PLOT FOLD COVERAGE COLUMN ----------------------------------------------------

cov_plt = fold_dat %>% 
  plot_sexchrom_scatter_two_tailed(xcol='mb', ycol='ratio', alpha=ALPHA) +
  coord_cartesian(y=c(-1,1)) +
  facet_wrap(~species, nrow=length(species_names))

# PLOT Y-LIKE DENSITY -----------------------------------------------------

#old one-tailed version
# ydens_plt = ydens_dat %>% 
#   plot_sexchrom_scatter(xcol='mb', ycol='ylike_density', ylim=c(0,0.0015), alpha=ALPHA) +
#   facet_wrap(~species, nrow=length(species_names))

ydens_plt = ydens_dat %>% 
  plot_sexchrom_scatter_two_tailed(xcol='mb', ycol='ylike_density', ylim=c(0,0.0015), alpha=ALPHA) +
  facet_wrap(~species, nrow=length(species_names))


# PLOT W-LIKE DENSITY -----------------------------------------------------

# #old one-tailed version
# wdens_plt = wdens_dat %>% 
#   mutate(species = factor(as.character(species),
#                           levels=c('wingei','picta', 'latipinna', 'Gambusia'))) %>% 
#   plot_sexchrom_scatter(xcol='mb', ycol='ylike_density', ylim=c(0,0.0015), alpha=ALPHA) +
#   facet_wrap(~species, nrow=4) +
#   labs(y='W-like density')


wdens_plt = wdens_dat %>% 
  mutate(species = factor(as.character(species),
                          levels=c('wingei','picta', 'latipinna', 'Gambusia'))) %>% 
  plot_sexchrom_scatter_two_tailed(xcol='mb', ycol='ylike_density', ylim=c(0,0.0015), alpha=ALPHA) +
  facet_wrap(~species, nrow=4) +
  labs(y='W-like density')



# DIFFERENTIAL EXPRESSION -------------------------------------------------

#Absolute log2 fold change
ge_plt=ge_res %>% 
  plot_sexchrom_scatter(xcol='mb', ycol='abslog2') +
  labs(y=bquote('absolute log'[2]~'M:F fold change')) +
  theme(axis.title.x = element_blank()) +
  facet_wrap(~species, nrow=length(species_names))

#differential expression
ge_plt = ge_res %>% 
  plot_sexchrom_scatter_two_tailed(xcol='mb', ycol='log2FoldChange', alpha=ALPHA) +
  labs(y=bquote('absolute log'[2]~'M:F fold change')) +
  facet_wrap(~species, nrow=length(species_names))


# VCF_WRAPPER COLUMN -------------------------------------------------------------
#didn't use these, but kept for reference

#ratio of male snp density to female (as in Darolti)
snp_dens_ratio = vcf_dat %>% 
  plot_sexchrom_scatter_two_tailed(xcol='mb', ycol='mf_snpCount_ratio',alpha=ALPHA) +
  labs(x='position (Mb)',
       y=bquote(log[2]~'M:F SNP density ratio')) +
  facet_wrap(~species, nrow=length(species_names)) 


#proportion male specific
vcf_dat %>% 
  plot_sexchrom_scatter_two_tailed(xcol='mb', ycol='pMaleSpecific', ylim=c(0,0.1)) +
  facet_wrap(~species, nrow=length(species_names)) +
  labs(y='proportion male-specific alleles')


#proportion female specific
female_specific = vcf_dat %>% 
  plot_sexchrom_scatter_two_tailed(xcol='mb', ycol='pFemaleSpecific', ylim=c(0,0.1), alpha=ALPHA) +
  labs(x='position (Mb)',
       y='proportion female-specific alleles') +
  facet_wrap(~species, nrow=length(species_names)) 


#mean FST
fst_plt = vcf_dat %>% 
  plot_sexchrom_scatter_two_tailed(xcol='mb', ycol='MEAN_FST', alpha=ALPHA) +
  labs(x='position (Mb)',
       y=bquote(mean~F[ST])) +
  facet_wrap(~species, nrow=length(species_names))


#Smale:female PI ratio
vcf_dat %>% 
  plot_sexchrom_scatter(xcol='mb', ycol='mf_pi_ratio') +
  facet_wrap(~species, nrow=length(species_names))
  
# ASSEMBLE COPY OF FIGURE 1 -----------------------------------------------

#choose the plots to build multipanel
plt_list = list(cov_plt,
                snp_dens_ratio,
                ydens_plt,
                female_specific,
                fst_plt,
                ge_plt)

title_list = list('Read\ndepth\nratio',
               'SNP\ndensity\nratio',
                    'SDR-like\nSNP\ndensity',
                    'Female-\nspecific\nallele density',
                    '\nFST\n',
                    'Expression\nratio\n')

#make modifications to the list for final plot
mod_plt_list = list()
for (i in 1:length(plt_list)){
  t = title_list[[i]]
  p = plt_list[[i]]
  mp = p + labs(title = t) +
    theme(axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5,
                              vjust = 0.5,
                              face = 'plain'))
  mod_plt_list[[i]] = mp
}
pans = plot_grid(plotlist = mod_plt_list, nrow=1)
xlab = ggdraw() + draw_label('position (Mb)')
plot_grid(pans, xlab, nrow=2, rel_heights = c(20,1))


# assign the SDR boundaries -----------------------------------------------

#select bound
SDR_start = 6.5
max_chr8 = ydens_dat %>% 
  filter(chr==8) %>% 
  pull(mb) %>% 
  max()
rev_SDR_start = max_chr8 - SDR_start

#rev
ydens_dat %>% 
  plot_sexchrom_scatter_two_tailed(xcol='mb', ycol='ylike_density', ylim=c(0,0.0015), alpha=ALPHA) +
  facet_wrap(~species, nrow=length(species_names)) +
  geom_vline(xintercept = rev_SDR_start)

#rev
fold_dat %>% 
  plot_sexchrom_scatter_two_tailed(xcol='mb', ycol='ratio', alpha=ALPHA) +
  coord_cartesian(y=c(-1,1)) +
  facet_wrap(~species, nrow=length(species_names)) +
  geom_vline(xintercept = rev_SDR_start)

save(SDR_start, rev_SDR_start, file='figure_plotting/sdr_start.Rdata')
  

#bonus things:
# assemble wingei for closer look -----------------------------------------

#plot ylike
yw = ydens_dat %>% 
  filter(species == 'wingei') %>% 
  mutate(species = factor(as.character(species),
                          levels=c('wingei'))) %>% 
  plot_sexchrom_scatter(xcol='mb', ycol='ylike_density', ylim=c(0,0.0015), alpha=ALPHA) +
  facet_wrap(~species, nrow=4) +
  labs(y='Y-like density')

#plot w-like
ww = wdens_dat %>% 
  filter(species == 'wingei') %>% 
  mutate(species = factor(as.character(species),
                          levels=c('wingei'))) %>% 
  plot_sexchrom_scatter(xcol='mb', ycol='ylike_density', ylim=c(0,0.0015), alpha=ALPHA) +
  facet_wrap(~species, nrow=4) +
  labs(y='W-like density')

#plot male-specific snps
msw = vcf_dat %>% 
  filter(species == 'wingei') %>% 
  mutate(species = factor(as.character(species),
                          levels=c('wingei'))) %>% 
  plot_sexchrom_scatter(xcol='mb', ycol='pMaleSpecific') +
  labs(x='position (Mb)',
       y='male-specific SNPs') +
  facet_wrap(~species, nrow=length(species_names)) 

#plot female-specific snps
fsw = vcf_dat %>% 
  filter(species == 'wingei') %>% 
  mutate(species = factor(as.character(species),
                          levels=c('wingei'))) %>% 
  plot_sexchrom_scatter(xcol='mb', ycol='pFemaleSpecific') +
  labs(x='position (Mb)',
       y='female-specific SNPs') +
  facet_wrap(~species, nrow=length(species_names)) 

#plot ratio
ratiow = vcf_dat %>% 
  filter(species == 'wingei') %>% 
  mutate(species = factor(as.character(species),
                          levels=c('wingei'))) %>% 
  plot_sexchrom_scatter(xcol='mb', ycol='mf_snpCount_ratio') +
  labs(x='position (Mb)',
       y=bquote('M:F SNP density')) +
  facet_wrap(~species, nrow=length(species_names)) 


#assemble
plot_grid(msw,
          fsw,
          yw,
          ww,
          ratiow,
          nrow=5,
          align='v')
  

