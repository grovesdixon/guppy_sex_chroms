
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

select_spp = 'picta'
ddat_in = paste('both_dna_rna/callXY/', select_spp, '_ylikeDepths.tsv', sep='')
ydat_in = paste('both_dna_rna/callXY/', select_spp, '_hcut8_y.vcf', sep='')
ddat = read_tsv(ddat_in)

# #remove failed samples
# ddat = ddat %>% 
#   select(-contains('picta_4Female_RNA'))

# normalize by total allele depths -----------------------------------------------

#get the names as in the vcf
vcf_samples0 = ddat %>% 
  select(contains('ref_')) %>% 
  colnames()
vcf_samples = sub('ref_', '', vcf_samples0)

#get totals by sample
depths = ddat %>% 
  select(contains(select_spp)) %>% 
  pivot_longer(everything()) %>% 
  mutate(name = sub('ref_', '', name),
         name = sub('alt_', '', name)) %>% 
  group_by(name) %>% 
  summarize(tot = sum(value))
tots = depths$tot
names(tots) = depths$name
refs = tots
names(refs) = paste('ref_', names(refs), sep='')
alts = tots
names(alts) = paste('alt_', names(alts), sep='')
norm_vec = append(refs, alts)

#normalize 
pos = ddat %>% 
  select(CHROM:ALT)
donly = data.frame(ddat)[,5:ncol(ddat)]
colnames(donly) == names(norm_vec)
norm_d = sweep(donly, 2, norm_vec/1e6, `/`)
norm_ddat = cbind(pos, norm_d) %>% 
  as_tibble

# merge with y-allele calls -----------------------------------------------

ydat = read_tsv(ydat_in) %>% 
  select(1:5) %>% 
  set_names(c('CHROM', 'POS', 'ID', 'Xallele', 'Yallele')) %>% 
  select(-ID)

#merge them
yddat = norm_ddat %>% 
  inner_join(ydat, by = c('CHROM', 'POS'))

#check we have all the sites as expected
nrow(yddat)==nrow(ddat)


# split by type -----------------------------------------------------------

#male dna
mdna = yddat %>% 
  select(contains('_DNA')) %>% 
  select(-contains('Female'))
#female dna
fdna = yddat %>% 
  select(contains('_DNA')) %>% 
  select(contains('Female'))
#male rna
mrna = yddat %>% 
  select(contains('_RNA')) %>% 
  select(-contains('Female'))
#female rna
frna = yddat %>% 
  select(contains('_RNA')) %>% 
  select(contains('Female'))
#xy calls
pos = yddat %>% 
  select(CHROM, POS, REF, ALT, Xallele, Yallele)
#assemble into list
datlist = list(mdna,
                fdna,
                mrna,
                frna)
names(datlist) = c('mdna',
                       'fdna',
                       'mrna',
                       'frna')

#get mean ref and alt for each
get_mean_by_starts_with = function(dat, start_string){
  sub = dat %>% 
    select(starts_with(start_string))
  mns = apply(sub, 1, function(x) mean(x, na.rm=TRUE))
  return(mns)
}

#function to get use the xy calls and the mean reference and 
#alternative allele depths to get mean X and Y depths
call_mn_xy = function(dat, n){
  pos %>% 
    mutate(mnRef = get_mean_by_starts_with(dat, 'ref_'),
           mnAlt = get_mean_by_starts_with(dat, 'alt_'),
           mnXallele = if_else(REF==Xallele,
                               mnRef,
                               mnAlt),
           mnYallele = if_else(REF==Yallele,
                               mnRef,
                               mnAlt)) %>% 
    select(-mnRef, -mnAlt) %>% 
    set_names('CHROM','POS','REF','ALT','Xallele','Yallele',
              paste(n, 'X', sep=''),
              paste(n, 'Y', sep=''))
}

#run call_mn_xy on each dataset
mn_datlist = list()
for (n in names(datlist)){
  mn_datlist[[n]] = call_mn_xy(datlist[[n]], n)
}

#assemble
ddat2 = purrr::reduce(mn_datlist, inner_join, by=c('CHROM','POS','REF','ALT','Xallele','Yallele')) %>% 
  mutate(sex_chrom = if_else(CHROM==8,
                             'sex',
                             'auto'))


save(ddat2, file=paste('both_dna_rna/callXY/', select_spp, '_ddatByDepth.Rdata', sep=''))

