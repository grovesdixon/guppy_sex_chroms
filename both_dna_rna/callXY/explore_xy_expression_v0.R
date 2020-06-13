#explore_xy_expression.R

rm(list=ls())
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

select_spp = 'picta'
select_spp = 'wingei'
select_spp = 'reticulata'


infile = paste('both_dna_rna/callXY/', select_spp, '_ddat.Rdata', sep='')
infile
ll=load(infile)



# total male vs female ------------------------------------------------------------

#for dna
tddat = ddat2 %>% 
  mutate(tm = mdnaX + mdnaY,
         tf = fdnaX + fdnaY)

#check accross chroms
tddat %>% 
  ggplot(aes(x=POS, y=log(tm/tf, 2))) +
  geom_hline(yintercept = 0, lty=1) +
  geom_point() +
  labs(y=bquote(log[2]*"M:F total SNP depth")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  facet_grid(~CHROM)

#plot against one another
tddat %>% 
  filter(CHROM==8) %>% 
  ggplot(aes(x=POS, y=log(tm/tf, 2))) +
  geom_hline(yintercept = 0, lty=1) +
  geom_point() +
  geom_smooth(se=FALSE) +
  labs(y=bquote(log[2]*"M:F total SNP depth"))

#plot against one another
tddat %>% 
  filter(CHROM==8) %>% 
  ggplot(aes(x=log(tm), y=log(tf))) +
  geom_point() +
  geom_abline(slope=1,intercept=0)



# total rna  --------------------------------------------------------------

#for dna
trdat = ddat2 %>% 
  mutate(trm = mrnaX + mrnaY,
         trf = frnaX + frnaY)

#check accross chroms
trdat %>% 
  ggplot(aes(x=POS, y=log(trm/trf, 2))) +
  geom_hline(yintercept = 0, lty=2) +
  geom_point() +
  # geom_smooth() +
  labs(y=bquote(log[2]*"M:F total SNP depth")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  facet_grid(~CHROM)


#plot against one another
trdat %>% 
  filter(CHROM==8) %>% 
  ggplot(aes(x=POS, y=log(trm/trf, 2))) +
  geom_hline(yintercept = 0, lty=2) +
  geom_point() +
  geom_smooth(se=FALSE) +
  labs(y=bquote(log[2]*"M:F total SNP depth"))


# total dna vs total rna --------------------------------------------------


tdat = trdat %>% 
  select(CHROM, POS, trm, trf) %>% 
  left_join(tddat) %>% 
  mutate(rna_ratio = log(trm/trf,2),
         dna_ratio = log(tm/tf,2))
  
#plot against each other
tdat %>% 
  ggplot(aes(x=dna_ratio,
             y=rna_ratio,
             color=sex_chrom)) +
  geom_abline(slope=1,intercept = 0) +
  geom_point() +
  scale_color_manual(values=c('black', 'red'))


#plot ratio of ratios
tdat %>% 
  mutate(diff = rna_ratio - dna_ratio) %>% 
  ggplot(aes(x=diff, fill=sex_chrom)) +
  geom_vline(xintercept = 0, lty=2) +
  geom_density(alpha=0.75) +
  labs(x='M:F RNA / M:F DNA')



# X vs Y stuff ------------------------------------------------------------

#male Y-like vs female Y-like 
#(females shoudl theoretically have 0 depth)
ddat2 %>% 
  select(mdnaY, fdnaY) %>% 
  pivot_longer(everything()) %>% 
  mutate(value=log(value, 2)) %>% 
  ggplot(aes(x=value, fill=name)) +
  geom_density(alpha=0.75) 



#X vs Y in males
malexy = ddat2 %>% 
  select(CHROM,POS,sex_chrom, starts_with('m')) %>% 
  mutate(dna_ratio = log(mdnaY / mdnaX, 2),
         rna_ratio = log(mrnaY / mrnaX, 2))

#scatterplot
malexy %>% 
  ggplot(aes(x=dna_ratio, y=rna_ratio, color=sex_chrom)) +
  geom_point()

#scatterplot
malexy %>% 
  filter(CHROM==8) %>% 
  ggplot(aes(x=dna_ratio, y=rna_ratio, fill=sex_chrom)) +
  geom_point()

#density for sex chromosome
malexy %>% 
  filter(CHROM==8) %>% 
  select(dna_ratio, rna_ratio) %>% 
  pivot_longer(everything()) %>% 
  mutate(value = if_else(value==-Inf,
                         -10,
                         value),
         value = if_else(value==Inf,
                         10,
                         value)) %>% 
  ggplot(aes(x=value, fill=name)) +
  geom_vline(xintercept = 0, lty=2) +
  geom_density(alpha=0.75)

#density for rest of genome
malexy %>% 
  filter(CHROM!=8) %>% 
  select(dna_ratio, rna_ratio) %>% 
  pivot_longer(everything()) %>% 
  mutate(value = if_else(value==-Inf,
                         -10,
                         value),
         value = if_else(value==Inf,
                         10,
                         value)) %>% 
  ggplot(aes(x=value, fill=name)) +
  geom_vline(xintercept = 0, lty=2) +
  geom_density(alpha=0.75)




# compensation barplot ----------------------------------------------------

x_dat = ddat2 %>% 
  mutate(DNA = log(mdnaX / fdnaX, 2),
         RNA = log(mrnaX / frnaX, 2)) %>% 
  select(CHROM, POS, sex_chrom, DNA, RNA)


x_dat %>% 
  pivot_longer(DNA:RNA,
               names_to = 'source',
               values_to = 'ratio') %>% 
  ggplot(aes(x=sex_chrom, y=ratio, fill=source)) +
  geom_hline(yintercept = 0, lty=2) +
  # geom_hline(yintercept = 1) +
  geom_boxplot() +
  labs(y='M:F X-linked allele depth')



# X depth over total depth ------------------------------------------------

x_depth = ddat2 %>% 
  mutate(male.DNA = mdnaX / (mdnaX + mdnaY),
         female.DNA = fdnaX / (fdnaX + fdnaY),
         male.RNA = mrnaX / (mrnaX + mrnaY),
         female.RNA = frnaX / (frnaX + frnaY)) %>% 
  select(CHROM, POS, sex_chrom, contains('male')) %>% 
  pivot_longer(contains('male')) %>% 
  separate(name, into=c('sex', 'source'))

#get stats
x_depth %>% 
  group_by(CHROM, sex, source) %>% 
  summarize(N=n()) %>% 
  data.frame()


#plot the proportion of depth from X chromosome
x_depth %>% 
  ggplot(aes(x=sex, y=value, fill=source)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.5, lty=2) +
  # geom_hline(yintercept = c(-1,1)) +
  # coord_cartesian(y=c(-2,2)) +
  labs(y='proportional X-like depth') +
  facet_grid(~sex_chrom)


#plot the same in terms of Y
x_depth %>% 
  mutate(yprop = 1-value) %>% 
  ggplot(aes(x=sex, y=yprop, fill=source)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty=2) +
  geom_hline(yintercept = 0.5, lty=2) +
  # coord_cartesian(y=c(-2,2)) +
  labs(y='proportional Y-like depth') +
  facet_grid(~sex_chrom)


#plot along chromosome
x_depth %>% 
  filter(CHROM==8,
         sex=='male') %>% 
  ggplot(aes(x=POS/1e6, y=value, color=source)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  labs(x='sex chromosome position (Mb)',
       y='proporitional X-like depth')

# stringent X depth over total depth ------------------------------------------------
#the Y proportion should be zero for females
ddat2 %>% 
  select(fdnaY, frnaY) %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x=value, fill=name)) + 
  geom_density() +
  lims(x=c(0,0.01))

x_depth = ddat2 %>% 
  filter(fdnaY==0,
         frnaY==0) %>% 
  mutate(male.DNA = mdnaX / (mdnaX + mdnaY),
         female.DNA = fdnaX / (fdnaX + fdnaY),
         male.RNA = mrnaX / (mrnaX + mrnaY),
         female.RNA = frnaX / (frnaX + frnaY)) %>% 
  select(CHROM, POS, sex_chrom, contains('male')) %>% 
  pivot_longer(contains('male')) %>% 
  separate(name, into=c('sex', 'source'))

#get stats
x_depth %>% 
  group_by(CHROM, sex, source) %>% 
  summarize(N=n()) %>% 
  data.frame()
  
  

#plot the proportion of depth from X chromosome
x_depth %>% 
  ggplot(aes(x=sex, y=value, fill=source)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.5, lty=2) +
  # geom_hline(yintercept = c(-1,1)) +
  # coord_cartesian(y=c(-2,2)) +
  labs(y='proportional X-like depth') +
  facet_grid(~sex_chrom)



# repeat in list mode -----------------------------------------------------

# load --------------------------------------------------------------------
species = c('reticulata',
            'wingei',
            'picta')


load_ddat = function(select_spp){
  print(select_spp)
  infile = paste('both_dna_rna/callXY/', select_spp, '_ddat.Rdata', sep='')
  ll=load(infile)
  return(ddat2)
}

ddat_list = map(species, load_ddat)
names(ddat_list) = species



# total DNA depth ---------------------------------------------------------

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


get_male_xy = function(dat){
  malexy = dat %>% 
    filter(mdnaY > 0) %>% 
    select(CHROM,POS,sex_chrom, starts_with('m')) %>% 
    mutate(DNA = log(mdnaY / mdnaX, 2),
           RNA = log(mrnaY / mrnaX, 2)) 
  return(malexy)
}

malexy_list = map(dat_list, get_male_xy)

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
         title = n)
  yx_dens_plts[[n]] = plt
}

plot_grid(plotlist = yx_dens_plts,
          nrow=3)



# proportion of coverage from X -------------------------------------------

#decide on female Y depth cutoff
#(aparently mpileup allows some depth for the alternative allele without calling it)
reduce(dat_list, rbind) %>% 
  ggplot(aes(x=mdnaY)) +
  geom_density() +
  lims(x=c(0,0.5))
  

FEMALE_DEPTH_CUT = 0.05 #female Y depth must be less than or equal to this (set to 10 to make irrelevent)
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

x_depth_list = map(dat_list, get_x_depth)

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
    geom_boxplot() +
    geom_hline(yintercept = 0.5, lty=2) +
    # geom_hline(yintercept = c(-1,1)) +
    # coord_cartesian(y=c(-2,2)) +
    labs(y='proportional X-like depth',
         subtitle = n) +
    facet_grid(~sex_chrom)
  comp_plts[[n]] = plt
}
plot_grid(plotlist = comp_plts,
          nrow=3)



save()
