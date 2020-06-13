#fold_coverage_analysis.R
rm(list=ls())
source('functions.R')
species_names

# SELECT FOLD COVERAGE DATASET ---------------------------------------------

#full maculatus
fc_folder = 'dna/fold_coverage/full_maculatus/'
infile = paste(fc_folder, 'all_multicov_results.tsv', sep='')
dat = read_tsv(infile, trim_ws = TRUE)
colnames(dat) = sub('Female', 'F_', colnames(dat))
colnames(dat) = sub('Male', 'M_', colnames(dat))
male_string = 'M_'
female_string = 'F_'


#quick maculatus
fc_folder = 'dna/fold_coverage/quick_maculatus/'
infile = paste(fc_folder, 'all_multicov_results.tsv', sep='')
dat = read_tsv(infile, trim_ws = TRUE)
male_string = 'M_'
female_string = 'F_'


#quick hellerii
fc_folder = 'dna/fold_coverage/quick_hellerii/'
infile = paste(fc_folder, 'all_multicov_results.tsv', sep='')
dat = read_tsv(infile, trim_ws = TRUE) %>% 
  swap_chrs()
male_string = 'M_'
female_string = 'F_'


#full hellerii
fc_folder = 'dna/fold_coverage/full_hellerii/'
infile = paste(fc_folder, 'all_multicov_results.tsv', sep='')
dat = read_tsv(infile, trim_ws = TRUE) %>% 
  swap_chrs()
male_string = 'M_'
female_string = 'F_'



# FORMAT ------------------------------------------------------------------

pos_dat = dat %>% 
  dplyr::select(chr:name)

cov_dat = dat %>% 
  select(-chr, -start, -end, -name)


# RUN FOR EACH IN A LOOP --------------------------------------------------

#set up identifying strings for the species in the column names
species_strings = c('reticulata',
                    'wingei',
                    'picta',
                    'latipinna',
                    'Gambusia')

fold_dat = data.frame()
for (i in 1:length(species_strings)){
  selection = species_strings[i]
  spp_name = species_names[i]
  print('-----------')
  print(selection)
  print(spp_name)
  
  #make the subset
  sdat = cov_dat %>% 
    dplyr::select(contains(selection))
  
  #get counts per million accross windows
  tots = apply(sdat, 2, function(x) sum(x, na.rm=TRUE))
  m = tots/1e6
  cpm_dat = sweep(sdat, 2, m, `/`) %>% 
    as_tibble()
  cpm_male = cpm_dat %>% 
    dplyr::select(contains(male_string))
  print('Male data:')
  print(head(cpm_male))
  print(dim(cpm_male))
  cpm_female = cpm_dat %>% 
    dplyr::select(contains(female_string))
  print('Female data:')
  print(head(cpm_female))
  print(dim(cpm_female))
  mn_male = apply(cpm_male, 1, function(x) mean(x, na.rm=TRUE))
  mn_female = apply(cpm_female, 1, function(x) mean(x, na.rm=TRUE))
  mn_all = apply(cpm_dat, 1, function(x) mean(x, na.rm=TRUE))
  mn_dat = pos_dat
  mn_dat$male = mn_male
  mn_dat$female = mn_female
  mn_dat$ratio = log((mn_male/mn_female), 2)
  mn_dat$all = mn_all
  #combine into single dataframe
  mn_dat$species = spp_name
  fold_dat = rbind(fold_dat, mn_dat)
}
head(fold_dat)
dim(fold_dat)
unique(fold_dat$species)
fold_dat = fold_dat %>% 
  mutate(species = as.character(species),
         mb = start/1e6) %>% 
  as_tibble()
outfile = paste(fc_folder, 'fold_dat.Rdata', sep='')
save(fold_dat, file=outfile)




# PLOT --------------------------------------------------------------------

cov_plt = fold_dat %>% 
  plot_sexchrom_scatter(xcol='mb', ycol='ratio', alpha=1) +
  facet_wrap(~species, nrow=length(species_names))
cov_plt



# PLOT EACH OF THEM -------------------------------------------------------

#load each of the four datasets
path_list = c('dna/fold_coverage/full_maculatus/fold_dat.Rdata',
              'dna/fold_coverage/quick_maculatus/fold_dat.Rdata',
              'dna/fold_coverage/full_hellerii/fold_dat.Rdata',
              'dna/fold_coverage/quick_hellerii/fold_dat.Rdata')
names(path_list) = c('full maculatus', '1 million maculatus', 'full hellerii', '1 million hellerii')
load_fd = function(path){
  load(path)
  return(fold_dat)
}
fd_list = map(path_list, load_fd)

#plot them
fd_plt_list = list()
for (n in names(fd_list)){
  fold_dat = fd_list[[n]] %>% 
    mutate(species = factor(species, levels = species_names))
  plt = fold_dat %>% 
    plot_sexchrom_scatter_two_tailed(xcol='mb', ycol='ratio', alpha=0.2) +
    labs(title = n,
         x = 'position (Mb)',
         y = bquote(log[2]~"M:F fold coverage")) +
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    facet_wrap(~species, nrow=length(species_names))
  fd_plt_list[[n]] = plt
}

# plot_grid(plotlist = fd_plt_list, nrow=1)

#just full
full= list(fd_plt_list[[1]] + labs(title='maculatus') + theme(plot.title = element_text(face='italic')),
           fd_plt_list[[3]] + labs(title = 'hellerii')+ theme(plot.title = element_text(face='italic')))
plot_grid(plotlist = full, nrow=1)


# MAKE SELECTION ----------------------------------------------------------


#reticulata
selection = 'reticulata'

#wingei
selection = 'P.win'

#picata
selection = 'P.pic'

#latipinna
selection = 'P.lat'

#G. holbrooki
selection = 'G.hol.'


#make the subset
sdat = cov_dat %>% 
  dplyr::select(contains(selection))

#get counts per million accross windows
tots = apply(sdat, 2, function(x) sum(x, na.rm=TRUE))
m = tots/1e6
cpm_dat = sweep(sdat, 2, m, `/`) %>% 
  as_tibble()
cpm_male = cpm_dat %>% 
  dplyr::select(contains(male_string))
cpm_female = cpm_dat %>% 
  dplyr::select(contains(female_string))
mn_male = apply(cpm_male, 1, function(x) mean(x, na.rm=TRUE))
mn_female = apply(cpm_female, 1, function(x) mean(x, na.rm=TRUE))
mn_all = apply(cpm_dat, 1, function(x) mean(x, na.rm=TRUE))
mn_dat = pos_dat
mn_dat$male = mn_male
mn_dat$female = mn_female
mn_dat$ratio = log((mn_male/mn_female), 2)
mn_dat$all = mn_all


ycol='ratio'
mn_dat %>% 
  mutate(mb=start/1e6) %>% 
  ggplot(aes_string(x='mb', y=ycol)) +
  geom_hline(yintercept = -1, color='red', lwd=0.5) +
  geom_point(alpha=0.2) +
  geom_smooth(lwd=0.5) +
  geom_hline(yintercept = 0, lty=2) +
  lims(y=c(-1.5, 1.5)) +
  labs(x='position (Mb)', y=bquote(log[2]~"M:F CPM")) +
  facet_wrap(~chr, scales = 'free_x')


#sex chrom
mn_dat %>% 
  filter(chr==8) %>% 
  mutate(mb=start/1e6) %>% 
  ggplot(aes_string(x='mb', y=ycol)) +
  geom_hline(yintercept = -1, color='red', lwd=0.5) +
  geom_point(alpha=0.4) +
  geom_smooth(lwd=0.5) +
  geom_hline(yintercept = 0, lty=2) +
  lims(y=c(-1.5, 1.5)) +
  labs(x='chromosome 8 (Mb)', y=bquote(log[2]~"M:F CPM"))


#get in long format
lcpm_dat = pos_dat %>% 
  dplyr::select(chr, start) %>% 
  cbind(cpm_dat) %>% 
  as_tibble() %>% 
  pivot_longer(contains(selection),
               names_to = 'sampleName',
               values_to = 'value') %>% 
  mutate(sex = if_else(grepl(male_string, sampleName),
                       'male',
                       'female'))

#plot lines by sex
lcpm_dat %>% 
  filter(chr==8) %>% 
  ggplot(aes(x=start, y=value)) +
  geom_smooth(aes(color=sampleName, lty=sex), se=FALSE)


#plot each individual compared to mean
head(mn_dat)
head(lcpm_dat)
mn_dat %>% 
  left_join(lcpm_dat, by = c('chr', 'start')) %>% 
  filter(chr==8) %>% 
  mutate(ratioToMean = log(value / all), 2,
         mb = start/1e6) %>% 
  ggplot(aes(x=mb, y=ratioToMean, color=sampleName, lty=sex)) +
  geom_hline(yintercept = 0, lty=2) +
  geom_smooth(se=FALSE) +
  labs(x='chr8 position (Mb)',
       y = bquote(log[2]~'ratio to mean CPM'))
