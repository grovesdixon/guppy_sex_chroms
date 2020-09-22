#explore_xy_expression.R
#for some messier, more in depth explorations see explore_xy_expression_v0.R

rm(list=ls())
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

# load --------------------------------------------------------------------
species = c('reticulata',
            'wingei',
            'picta')


#function to load results from the format_spp_byTotalDepth.R scripts
load_ddat = function(select_spp){
  print(select_spp)
  infile = paste('both_dna_rna/callXY/', select_spp, '_ddatByDepth.Rdata', sep='')
  ll=load(infile)
  return(ddat2)
}

#apply the function
ddat_list = map(species, load_ddat)
names(ddat_list) = species
map(ddat_list, head)
names(ddat_list)

#now we have the SDR-like SNPs in a dataframe for each species

# total DNA depth ---------------------------------------------------------
#here we are looking at the total depth based on the summed allele depths
#this is largely to confirm the lower DNA depth results seen previously for P. picta males indicated degeration of Y

#check accross chroms
get_total_dna = function(dat){
  #for dna
  tddat = dat %>% 
    mutate(tm = mdnaX + mdnaY,
           tf = fdnaX + fdnaY)
  return(tddat)
}

#build plots
tddat_list = map(ddat_list, get_total_dna)
td_plt_list = list()
for (n in names(tddat_list)){
  tddat = tddat_list[[n]]
  plt = tddat %>% 
    ggplot(aes(x=POS, y=log(tm/tf, 2))) +
    geom_hline(yintercept = 0, lty=2) +
    geom_point(alpha = 0.1) +
    labs(y=bquote(log[2]*"M:F total SNP depth"),
         title = n) +
    lims(y=c(-2,2)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank()) +
    facet_grid(~CHROM)
  td_plt_list[[n]] = plt
}

#plot
plot_grid(plotlist = td_plt_list,
          nrow=3)

#repeat for sex chrom only
td_plt_list = list()
for (n in names(tddat_list)){
  tddat = tddat_list[[n]] %>% 
    filter(CHROM==8)
  plt = tddat %>% 
    ggplot(aes(x=POS, y=log(tm/tf, 2))) +
    geom_hline(yintercept = 0, lty=2) +
    geom_point(alpha = 0.1) +
    labs(y=bquote(log[2]*"M:F total SNP depth"),
         title = n) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank())
  td_plt_list[[n]] = plt
}

#plot
plot_grid(plotlist = td_plt_list,
          nrow=3)



# Y to X ratio in males ---------------------------------------------------
#now get the ratio of Y-allele depth and X-allele depth in males

get_male_xy = function(dat){
  malexy = dat %>% 
    filter(mdnaY > 0) %>% 
    select(CHROM,POS,sex_chrom, starts_with('m')) %>% 
    mutate(DNA = log(mdnaY / mdnaX, 2),
           RNA = log(mrnaY / mrnaX, 2)) 
  return(malexy)
}

malexy_list = map(ddat_list, get_male_xy)

#density for sex chromosome
yx_dens_plts = list()
for (n in names(malexy_list)){
  malexy = malexy_list[[n]]
  plt=malexy %>% 
    filter(CHROM==8) %>% 
    select(DNA, RNA) %>% 
    pivot_longer(everything()) %>% 
    mutate(value = if_else(value==-Inf,
                           -10,
                           value),
           value = if_else(value==Inf,
                           10,
                           value)) %>% 
    ggplot(aes(x=value, fill=name)) +
    geom_vline(xintercept = 0, lty=2) +
    geom_density(alpha=0.75) +
    labs(fill='',
         x = bquote(log[2]*'Y:X allele depth ratio'),
         subtitle = n) +
    theme(plot.subtitle = element_text(face='italic'))
  yx_dens_plts[[n]] = plt
}

plot_grid(plotlist = yx_dens_plts,
          nrow=3)



# SNPs with zero coverage for Ys from DNA ---------------------------------
#this is strange. There are aparently about 100 loci where 
#males were genotyped as heterozygous with allele depths of zero.

w = ddat_list[[2]]
w %>% 
  filter(mdnaY == 0)



# proportion of coverage from X -------------------------------------------

#decide on female Y depth cutoff
#(aparently mpileup allows some depth for the alternative allele without calling it)
purrr::reduce(ddat_list, rbind) %>% 
  ggplot(aes(x=mdnaY)) +
  geom_density() +
  lims(x=c(0,0.5))
  

FEMALE_DEPTH_CUT = 0.05 #female Y depth must be less than or equal to this (set to 10 to make irrelevent)
FEMALE_DEPTH_CUT = 10 #female Y depth must be less than or equal to this (set to 10 to make irrelevent)
get_x_depth = function(dat){
  x_depth = dat %>% 
    filter(fdnaY <= FEMALE_DEPTH_CUT,
           frnaY <= FEMALE_DEPTH_CUT) %>% 
    mutate(male.DNA = mdnaX / (mdnaX + mdnaY),
           female.DNA = fdnaX / (fdnaX + fdnaY),
           male.RNA = mrnaX / (mrnaX + mrnaY),
           female.RNA = frnaX / (frnaX + frnaY)) %>% 
    select(CHROM, POS, sex_chrom, contains('male')) %>% 
    pivot_longer(contains('male')) %>% 
    separate(name, into=c('sex', 'source'))
  return(x_depth)
}

x_depth_list = map(ddat_list, get_x_depth)



#get stats
for (n in names(x_depth_list)){
  print('----------')
  print(n)
  x_depth = x_depth_list[[n]]
  r=x_depth %>% 
    group_by(CHROM, sex, source) %>% 
    summarize(N=n()) %>% 
    data.frame()
  print(r)
}

#plot the proportion of depth from X chromosome
comp_plts = list()
for (n in names(x_depth_list)){
  x_depth = x_depth_list[[n]]
  plt = x_depth %>% 
    ggplot(aes(x=sex, y=value, fill=source)) +
    geom_hline(yintercept = 0.5, lty=2) +
    geom_boxplot() +
    # geom_hline(yintercept = c(-1,1)) +
    # coord_cartesian(y=c(-2,2)) +
    labs(y='proportional X-like depth',
         subtitle = n) +
    theme(plot.subtitle = element_text(face = 'italic'),
          axis.title = element_blank()) +
    facet_grid(~sex_chrom)
  comp_plts[[n]] = plt
}

l = cowplot::get_legend(comp_plts[[1]] + 
                          theme(legend.position = 'bottom',
                                legend.justification = 0.5))
comp_plts = lapply(comp_plts, function(x) return(x + theme(legend.position='none')))
pans = plot_grid(plotlist = comp_plts,
                 nrow=3)
ylab = ggdraw() + draw_label('proportion of coverage from X', angle=90)
top = plot_grid(ylab, pans, nrow=1, rel_widths = c(1,15))
plot_grid(top, l, nrow=2, rel_heights = c(25,1))

