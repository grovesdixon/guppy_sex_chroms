#compare_vcf_depths_with_featurecounts.R
#idea here is to ensure that a gene's mean depth from a vcf matches
#reasonably well with its expression estimate based on featurecounts.



library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())



#read in the vcf mean depths
vdat = read_tsv('both_dna_rna/vcf_depths/depthsByGene_picta.tsv')

#get totals
rvdat = vdat %>% 
  select(starts_with('ref'))
avdat = vdat %>% 
  select(starts_with('alt'))
colnames(rvdat) = sub('ref_', '', colnames(rvdat))
colnames(avdat) = sub('alt_', '', colnames(avdat))
sum(colnames(rvdat)==colnames(avdat)) == ncol(rvdat)
tvdat = rvdat + avdat
tvdat$Geneid = vdat$gene

#read in the featurecounts
ll=load('both_dna_rna/full_featurecounts.Rdata')
ll

#subset for same
colnames(tvdat) %in% colnames(full_dat)
matched = full_dat %>% 
  select(colnames(tvdat))


#convert to cpm
get_cpm = function(dat){
  counts = dat %>% 
    select(-c('Geneid'))
  m = apply(counts, 2, function(x) sum(x, na.rm=TRUE)) / 1e6
  cpm = sweep(counts, 2, m, `/`)
  cpm$Geneid = dat$Geneid
  cpm %>% 
    as_tibble()
}


fcpm = get_cpm(matched)
vcpm =get_cpm(tvdat)

fcpm$method = 'featurecounts'
vcpm$method = 'vcf'

lfdat = fcpm %>% 
  pivot_longer(-c('Geneid', 'method'),
               names_to = 'sample',
               values_to = 'fcpm')

lvdat = vcpm %>% 
  pivot_longer(-c('Geneid', 'method'),
               names_to = 'sample',
               values_to = 'vcpm')


combined = lfdat %>% 
  left_join(lvdat, by = c('Geneid', 'sample'))


#plot log counts
combined %>% 
  # filter(!sample %in% c('picta_4Female_RNA', 'picta_6Male_RNA')) %>% 
  filter(vcpm > 0,
         fcpm > 0) %>% 
  ggplot(aes(x=log(fcpm, 2), y = log(vcpm, 2), color=sample)) +
  geom_point() +
  geom_smooth()


#plot mean log counts
mns = combined %>% 
  # filter(!sample %in% c('picta_4Female_RNA', 'picta_6Male_RNA')) %>% 
  filter(vcpm > 0,
         fcpm > 0) %>% 
  group_by(Geneid) %>% 
  summarize(mnf = log(mean(fcpm, na.rm=TRUE), 2),
            mnv = log(mean(vcpm, na.rm=TRUE), 2))
lm1 = lm(mns$mnv ~mns$mnf)
summary(lm1)
mns %>% 
  ggplot(aes(x=mnf, y = mnv)) +
  geom_point() +
  geom_smooth()


#not fantastic, but seems workable
