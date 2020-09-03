#y_like_density.R
rm(list=ls())
source('functions.R')
species_names



# select dataset from call_xy.R --------------------------------------------

#the W-like frequency
species_names = c("wingei","picta","latipinna","Gambusia")
inputs = paste('both_dna_rna/callXY/', species_names, '_hcut8_w.vcf', sep='') #the W-like frequency
output_name = 'wdens_dat.Rdata'


#or

#the Y-like frequency
species_names = c("reticulata", "wingei","picta","latipinna","Gambusia" )
inputs = paste('both_dna_rna/callXY/', species_names, '_hcut8_y.vcf', sep='') #the Y-like frequency
output_name = 'ydens_dat.Rdata'



# load the output from call_xy.R --------------------------------------------
#
load_xy_vcf = function(input_file){
  read_tsv(input_file) %>% 
    dplyr::select(`#CHROM`, POS) %>% 
    set_names(c('chr', 'pos'))
}

vdat_list = map(inputs, load_xy_vcf)
names(vdat_list) = species_names



# get density along windows -----------------------------------------------

get_ylike_density = function(dat, window_size, chr_lengths){
  chrs = names(chr_lengths)
  res = data.frame()
  for (ichr in chrs){
    csub = dat %>% 
      filter(chr==ichr)
    chr_len = chr_lengths[ichr]
    breaks = seq(1, chr_len, window_size)
    breaks[length(breaks)] = chr_len
    h = hist(csub$pos,
             breaks = breaks,
             plot = FALSE)
    sub_res = data.frame(chr = ichr,
                         starts = h$breaks[1:(length(h$breaks)-1)],
                         ends = h$breaks[2:length(h$breaks)],
                         nYlike = h$counts)
    res = rbind(res, sub_res)
  }
  return(res)
}



# compare the y-like density between chromes ------------------------------

#load chr lengths for making bins
ll = load('metadata/maculatus_chr_lengths.Rdata')
ll

#make df for merging
cdat = data.frame(chr_lengths) %>% 
  rownames_to_column('chr') %>% 
  mutate(chr = as.numeric(chr))

#build a barplot
total_dens_list = list()
for (ispp in names(vdat_list)){
  vdat = vdat_list[[ispp]]
  ddat = vdat %>% 
    group_by(chr) %>% 
    summarize(n_ylike = n()) %>% 
    left_join(cdat, by = 'chr') %>% 
    mutate(density = n_ylike / chr_lengths,
           species = ispp)
  total_dens_list[[ispp]] = ddat
}
total_dens_df = purrr::reduce(total_dens_list, rbind)

selection = c('reticulata', 'wingei', 'picta')
total_dens_df %>% 
  filter(species %in% selection) %>% 
  mutate(sex = if_else(chr == 8,
                       'sex',
                       'autosome'),
         species = factor(species, levels=selection),
         chr = factor(chr, levels = 1:24)) %>% 
  ggplot(aes(x = chr, y=density, fill = sex)) +
  geom_bar(stat='identity') +
    labs(fill='') +
  facet_wrap(~species, scales = 'free_y', nrow=3)



# build alternative to barplot that's just sex vs autosomes ---------------

#boxplots
total_dens_df$species[total_dens_df$species=='reticulata'] <- 'P. reticulata'
total_dens_df$species[total_dens_df$species=='wingei'] <- 'P. wingei'
total_dens_df$species[total_dens_df$species=='picta'] <- 'P. picta'
total_dens_df %>% 
  filter(species %in% c('P. reticulata', 'P. wingei', 'P. picta')) %>% 
  mutate(sex_chrom = if_else(chr==8,
         'sex',
         'auto'),
         species = factor(species, levels = c('P. reticulata', 'P. wingei', 'P. picta')),
         sex_chrom = factor(sex_chrom, levels=c('auto', 'sex'))) %>% 
  ggplot(aes(x=sex_chrom, y = density, color=sex_chrom)) +
  geom_boxplot() +
  labs(x='chromosome') +
  scale_color_manual(values=c('grey50', 'black')) +
  facet_wrap(~species, scales = 'free') +
  theme(strip.text = element_text(face = 'italic'))
  


  


# look at windows ---------------------------------------------------------

#set window size
window_size = 1e4 #10Kb
window_size = 1e5 #100Kb

#run for each species
ydens_list = list()
for (spp in species_names){
  dat = vdat_list[[spp]]
  spp_ydens = get_ylike_density(dat, window_size, chr_lengths) %>% 
    mutate(species = spp)
  ydens_list[[spp]] = spp_ydens
}

#merge into single df
ydens_dat = purrr::reduce(ydens_list, rbind) %>% 
  as_tibble

#check things make sens
ydens_dat %>% 
  group_by(species) %>% 
  summarize(totYlike = sum(nYlike))


w = ydens_dat %>% 
  filter(species=='wingei') 

hist(w$starts)



# save restuls ------------------------------------------------------------

out_path = paste('both_dna_rna/callXY/', output_name, sep='')
print(out_path)
save(ydens_dat, file=paste('both_dna_rna/callXY/', output_name, sep=''))

