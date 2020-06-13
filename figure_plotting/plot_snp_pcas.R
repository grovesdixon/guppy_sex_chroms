

#choose pca files to input
pattern = 'filt2.vcf'
file_list = list.files(path = 'rnaseq/pcas/maculatus',
                       pattern = pattern,
                       full.names = TRUE)
file_names = list.files(path = 'rnaseq/pcas/maculatus',
                        pattern = pattern)
names(file_list) = file_names

#function to load pca object output from basic_snp_pca.R
load_pca = function(file_path){
  print(paste('loading ', file_path, '...', sep=''))
  ll=load(file_path)
  pca_df = pca$scores %>% 
    data.frame() %>% 
    rownames_to_column('sample') %>% 
    as_tibble()
  return(pca_df)
}

#load
pca_datlist = map(file_list, load_pca)


#function to plot
plot_pca = function(pca_df, pc1, pc2){
  pca_df %>% 
    mutate(spp = substr(sample, start=1, stop = 5),
           spp = if_else(grepl('^SRR', spp), 'P.reticulata', spp),
           sex = if_else(grepl('F', sample),
                         'female',
                         'male')) %>% 
    ggplot(aes_string(x=pc1, y=pc2, color='spp', shape='sex')) +
    geom_point(size=4)
}

names(pca_datlist)
plt_list = list()
for (n in names(pca_datlist)){
  plt_list[[n]] = plot_pca(pca_datlist[[n]], 'PC1', 'PC2') + labs(subtitle=n)
  plt_list[[paste(n,2,sep='')]] = plot_pca(pca_datlist[[n]], 'PC1', 'PC3') + labs(subtitle=n)
}
l = get_legend(plt_list[[1]] + theme(legend.position = 'bottom'))
plt_list2 = map(plt_list, function(x) return(x+theme(legend.position = 'none')))
pans = plot_grid(plotlist = plt_list2, nrow=3)
plot_grid(pans, l, nrow=2, rel_heights = c(12,1))



# LOOK AT WEIRD PICATA SSAMPLE --------------------------------------------
pca_df = pca_datlist[['chr2_filt2.vcf_pca.Rdata']]
pca_df %>% 
  filter(grepl('P.pic', sample),
         PC3 < 20,
         PC3 > -20)

#this is the same one with the weird heterozygosity
hdat = read_tsv('rnaseq/observed_heterozygosities/picta_heterozygosities.tsv')
read_heterozygosity = function(file_path){
  hdat = read_tsv(file_path) %>% 
    group_by(INDV) %>% 
    summarize(mn_obs_het = 1-mean(`O(HOM)`/N_SITES))
}

read_heterozygosity('rnaseq/observed_heterozygosities/picta_heterozygosities.tsv') %>% 
  ggplot(aes(x=INDV, y=mn_obs_het)) +
  geom_bar(stat='identity') +
  labs(x='P. picta individual',
       y='Mean observed heterozygosity')

#compare with ALL heterozygosities
h_files = list.files(path='rnaseq/observed_heterozygosities/',
                     pattern = 'heterozygosities.tsv',
                     full.names = TRUE)
h_datlist = map(h_files, read_heterozygosity)
purrr::reduce(h_datlist, rbind) %>% 
  mutate(species = substr(INDV, start=1, stop=4)) %>% 
  ggplot(aes(x=INDV, y=mn_obs_het, fill = species)) +
  geom_bar(stat='identity') +
  labs(x='P. picta individual',
       y='Mean observed heterozygosity') +
  theme(axis.text.x = element_text(angle=45))

#look at heterozygosity from DNA
read_heterozygosity('dna/observed_heterozygosities/picta_heterozygosities.tsv') %>% 
  ggplot(aes(x=INDV, y=mn_obs_het)) +
  geom_bar(stat='identity') +
  labs(x='P. picta individual',
       y='Mean observed heterozygosity') +
  coord_cartesian(ylim = c(0,0.8))

#the difference isn't there!!!!!!!!!

#doublechech the catting commands from cat_guppy_trim_files.txt (made with overview_data.R)
# cat SRR8785851_1.trim SRR8785852_1.trim SRR8785853_1.trim SRR8785854_1.trim > P.pic._2F_T_1_RNA-Seq_catted.fq
# cat SRR8785851_2.trim SRR8785852_2.trim SRR8785853_2.trim SRR8785854_2.trim > P.pic._2F_T_2_RNA-Seq_catted.fq
#to be absolutely sure, I re-download the SRA table for PRJNA528814 and searched for the ID again.
#Results are pasted here:
# Assay_Type	Run	Sample_Name
# RNA-Seq	SRR8785851	P.pic._2F_T
# RNA-Seq	SRR8785852	P.pic._2F_T
# RNA-Seq	SRR8785853	P.pic._2F_T
# RNA-Seq	SRR8785854	P.pic._2F_T











