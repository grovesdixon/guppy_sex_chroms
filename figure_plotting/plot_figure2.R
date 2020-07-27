#plot_figure2.R
#plot gene expression figures, including X:Y comparisons and male:female

rm(list=ls())
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

# load --------------------------------------------------------------------
species = c('reticulata',
            'wingei',
            'picta')

#These input files are built with format_[speciesname]_byTotalDepth.R
load_ddat = function(select_spp){
  print(select_spp)
  infile = paste('both_dna_rna/callXY/', select_spp, '_ddatByDepth.Rdata', sep='')
  ll=load(infile)
  return(ddat2)
}

ddat_list = map(species, load_ddat)
names(ddat_list) = paste('P.', species)

# Y to X ratio in males ---------------------------------------------------


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
    theme(plot.subtitle = element_text(face='italic',
                                       hjust = 0.5),
          legend.position = 'none',
          axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  yx_dens_plts[[n]] = plt
}

#format the panels
yx_dens_plts[[1]] = yx_dens_plts[[1]] +
  theme(legend.position = c(0.7,0.8))
ylab = ggdraw() + draw_label('Density', angle = 90)
xlab =ggdraw() + draw_label(bquote(log[2]*'Y:X allele depth ratio'))
pans = plot_grid(plotlist = yx_dens_plts,
          nrow=1)
right = plot_grid(pans, xlab, nrow=2, rel_heights = c(12,1))
full = plot_grid(ylab, right, nrow=1, rel_widths = c(1,20))
full





# plot male vs female boxplots ---------------------------------------------

#load the data (build these with both_dna_rna/featureCounts/explore_featurecounts.R)
ll=load('figure_plotting/fc_boxplot_objects.Rdata')
ll

#function to plot boxplots
plot_fc_boxes = function(fc_dat, spp){
  fc_dat %>% 
    mutate(sex_chrom = if_else(Chr==8,
                               'sex',
                               'auto')) %>% 
    ggplot(aes(x=sex_chrom, y=log2FoldChange, fill=source)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, lty=2) +
    labs(y=bquote(log[2]~'M:F coverage'),
         subtitle = spp,
         fill='') +
    coord_cartesian(y=c(-2, 2)) +
    theme(plot.subtitle = element_text(face = 'italic',
                                       hjust=0.5),
          axis.title = element_blank(),
          legend.position = 'bottom')
}

#build the plots
l = cowplot::get_legend(pplt)
lplt = plot_grid(l)
rplt = plot_fc_boxes(ret_fcs2, 'P. reticulata')
wplt = plot_fc_boxes(wing_fcs2, 'P. wingei') + theme(axis.text.y = element_blank())
pplt = plot_fc_boxes(picta_fcs2, 'P. picta') + theme(axis.text.y = element_blank())
plt_list = list(rplt, wplt, pplt)
plt_list = lapply(plt_list, function(x) return(x+theme(legend.position = 'none')))
pans = plot_grid(plotlist = plt_list, nrow=1)

#asemble
ylab = ggdraw() + draw_label(bquote(log[2]~'M:F read depth'), angle=90)
top = plot_grid(ylab, pans, nrow=1, rel_widths = c(1,20))
full = plot_grid(lplt, top, nrow=2, rel_heights = c(1,20))
full




# breakdown picta and wingei by SDR ---------------------------------------

#load eyeballed start positions for SDR 
ll=load('figure_plotting/sdr_start.Rdata')
ll
SDR_start #in Mb
rev_SDR_start #in Mb (for when coodinates are reversed)

#function to plot boxplots for SDR and PAR
plot_fc_boxes_par_sdr = function(fc_dat, spp){
  fc_dat %>% 
    mutate(sex_chrom = if_else((Chr==8) & Start/1e6 < SDR_start,
                               'PAR',
                               'auto'),
           sex_chrom = if_else((Chr==8) & Start/1e6 > SDR_start,
                               'SDR',
                               sex_chrom)) %>% 
    ggplot(aes(x=sex_chrom, y=log2FoldChange, fill=source)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, lty=2) +
    labs(y=bquote(log[2]~'M:F coverage'),
         subtitle = spp,
         fill='') +
    coord_cartesian(y=c(-2, 2)) +
    theme(plot.subtitle = element_text(face = 'italic',
                                       hjust=0.5),
          axis.title = element_blank(),
          legend.position = 'bottom')
}

#build the plots
wplt = plot_fc_boxes_par_sdr(wing_fcs2, 'P. wingei')
pplt = plot_fc_boxes_par_sdr(picta_fcs2, 'P. picta') + theme(axis.text.y = element_blank())
l = cowplot::get_legend(wplt)
lplt = plot_grid(l)
plt_list = list(wplt, pplt)
plt_list = lapply(plt_list, function(x) return(x+theme(legend.position = 'none')))
pans = plot_grid(plotlist = plt_list, nrow=1)

#asemble
ylab = ggdraw() + draw_label(bquote(log[2]~'M:F read depth'), angle=90)
top = plot_grid(ylab, pans, nrow=1, rel_widths = c(1,20))
full = plot_grid(lplt, top, nrow=2, rel_heights = c(1,20))
full

