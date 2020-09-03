#explore_featurecounts.R
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
rm(list=ls())

#load the results from featureCounts (see GET COUNTS FOR BASIC GENE EXPRESSION)
fc_path = 'both_dna_rna/featureCounts/'
infiles = list.files(path = fc_path, pattern = '*.tsv')

raw_datlist = list()
for (i in infiles){
  chr = strsplit(i, '_')[[1]][1]
  chr_num = as.numeric(sub('chr', '', chr))
  subchr = paste(chr, '_', sep='')
  infile = paste(fc_path, i, sep='')
  dat = read_tsv(infile, comment='#') %>% 
    filter(Chr==chr_num)
  colnames(dat) = sub('.bam', '', colnames(dat), fixed=TRUE)
  colnames(dat) = sub(subchr, '', colnames(dat), fixed=TRUE)
  raw_datlist[[chr]] = dat
  print('-----')
  print(paste(i,'...',sep=''))
  print(chr)
  print(ncol(dat))
  print(unique(dat$Chr))
}
map(raw_datlist, ncol)
full_dat = purrr::reduce(raw_datlist, rbind) %>% 
  as_tibble
head(full_dat)
save(full_dat, file='both_dna_rna/full_featurecounts.Rdata')

#check totals
samples = colnames(full_dat)[7:ncol(full_dat)]
spp = sapply(samples, function(x) strsplit(x, '_')[[1]][1])
sex = grepl('Female', sapply(samples, function(x) strsplit(x, '_')[[1]][2]))
sex[sex == TRUE]<-'Female'
sex[sex==FALSE]<-'Male'
molec = grepl('_DNA', samples)
molec[molec == TRUE]<-'DNA'
molec[molec==FALSE]<-'RNA'
spp_sex = paste(spp, sex, sep = '_')
spp_sex_molec = paste(spp_sex, molec)
table(spp)
table(sex)
table(spp_sex)
table(spp_sex_molec)

# convert to cpm ---------------------------------------------------------

pos = full_dat %>% 
  dplyr::select(Geneid:Length)
counts = full_dat %>% 
  dplyr::select(7:ncol(full_dat))

#normalize basica way
no_sex_counts = counts[pos$Chr != 8,]
m = apply(no_sex_counts, 2, function(x) sum(x, na.rm=TRUE)/1e6)
cpm = sweep(counts, 2, m, `/`)

# #use deseq to normalize
library(DESeq2)
coldata = data.frame(colnames(counts))
dds<-DESeqDataSetFromMatrix(data.frame(counts),
                            colData = coldata,
                            design = formula(~1))
rld = vst(dds)
cpm = data.frame(assay(rld))

cpm_dat = cbind(pos, cpm) %>% 
  as_tibble

# split by species --------------------------------------------------------

species = c('reticulata',
            'wingei',
            'picta',
            'latipinna',
            'Gambusia')


fcdat_list = list()
for (spp in species){
  sub = cpm_dat %>% 
    dplyr::select(Geneid:Length, starts_with(spp))
  fcdat_list[[spp]] = sub
}

# split by assay ----------------------------------------------------------

dna_list = list()
rna_list = list()

for (spp in species){
  spdat = fcdat_list[[spp]]
  dna = spdat %>% 
    dplyr::select(Geneid:Length,ends_with('DNA'))
  rna = spdat %>% 
    dplyr::select(Geneid:Length,ends_with('RNA'))
  dna_list[[spp]] = dna
  rna_list[[spp]] = rna
}

#check
map(dna_list, colnames)
map(rna_list, colnames)

# compare male and female coverage ----------------------------------------

compare_male_female = function(dat){
  pos = dat %>% 
    dplyr::select(Geneid:Length)
  male = dat %>% 
    dplyr::select(contains('Male')) %>% 
    dplyr::select(-contains('Female'))
  female = dat %>% 
    dplyr::select(contains('Female_'))
  m_means = apply(male, 1, function(x) median(x, na.rm=TRUE))
  f_means = apply(female, 1, function(x) median(x, na.rm=TRUE))
  mf_ratio = m_means/f_means
  pos$mf_ratio = mf_ratio
  return(pos)
}

mf_dna = map(dna_list, compare_male_female)
mf_rna = map(rna_list, compare_male_female)


spp = 'picta'
sex_chrom = 8

d = mf_dna[[spp]] %>% 
  filter(Chr == sex_chrom)
r = mf_rna[[spp]] %>% 
  filter(Chr == sex_chrom)


d %>% 
  mutate(lmf = log(mf_ratio, 2)) %>% 
  ggplot(aes(x=Start, y = lmf)) +
  geom_point() +
  geom_smooth()

r %>% 
  mutate(lmf = log(mf_ratio, 2)) %>% 
  ggplot(aes(x=Start, y = lmf)) +
  geom_point() +
  geom_smooth()


#plot compensation barplot
fold_comp_plts = list()
for (n in c('reticulata',
            'wingei',
            'picta')){
  d = mf_dna[[n]] %>% 
    mutate(source = 'DNA')
  r = mf_rna[[n]] %>% 
    mutate(source='RNA')
  plt = rbind(d,r) %>% 
    mutate(sex_chrom = if_else(Chr==8,
                               'sex',
                               'autosomes')) %>% 
    ggplot(aes(x=sex_chrom, y=log(mf_ratio, 2), fill=source)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, lty=2) +
    # geom_hline(yintercept = c(-1,1)) +
    # coord_cartesian(y=c(-2,2)) +
    labs(y=bquote(log[2]~'M:F coverage'),
         subtitle = n) +
    coord_cartesian(y=c(-2, 2)) +
    theme(plot.subtitle = element_text(face = 'italic'),
          axis.title.x = element_blank())
  fold_comp_plts[[n]] = plt
}
l = cowplot::get_legend(fold_comp_plts[[1]] + 
                          theme(legend.position = 'bottom',
                                legend.justification = 0.5))
fold_comp_plts = lapply(fold_comp_plts, function(x) return(x + theme(legend.position='none',
                                                                     axis.title.y = element_blank())))
pans = plot_grid(plotlist = fold_comp_plts,
          nrow=3)
ylab = ggdraw() + draw_label(bquote(log[2]~'M:F coverage'), angle=90)
top = plot_grid(ylab, pans, nrow=1, rel_widths = c(1,15))
plot_grid(top, l, nrow=2, rel_heights = c(25,1))


# use deseq to make M:F comparisons ---------------------------------------
library(DESeq2)
sub_counts = function(spp, assay){
  full_dat %>% 
    select(Geneid, contains(spp)) %>% 
    select(Geneid, contains(assay)) %>% 
    column_to_rownames('Geneid')
}
build_coldata = function(dat){
  samples = colnames(dat)
  coldata = data.frame(samples,
                       sex = if_else(grepl('Female',samples),
                                     'female',
                                     'male'))
  return(coldata)
}
sex_genes = pos %>% 
  filter(Chr == 8) %>% 
  pull(Geneid)
run_mf_deseq = function(count_dat,
                        coldat){
  no_sex = count_dat[!rownames(count_dat) %in% sex_genes,]
  dds_nosex<-DESeqDataSetFromMatrix(no_sex,
                                    colData = coldat,
                                    design = ~sex)
  dds<-DESeqDataSetFromMatrix(count_dat,
                              colData = coldat,
                              design = ~sex)
  dds_nosex<-estimateSizeFactors(dds_nosex)
  sizeFactors(dds) <- sizeFactors(dds_nosex)
  dds<-estimateDispersions(dds, fitType = 'local')
  dds<-nbinomWaldTest(dds)
  res = results(dds, contrast = c('sex', 'male', 'female'))
  return(res)
}

format_res = function(res,
                      spp,
                      assay){
  res %>% 
    data.frame() %>% 
    rownames_to_column('Geneid') %>% 
    left_join(pos, by = 'Geneid') %>% 
    mutate(species = spp,
           source = assay) %>% 
    as_tibble()
}

#RUN FOR PICTA
picta_dna = sub_counts('picta', 'DNA')
picta_rna = sub_counts('picta', 'RNA')
picta_dna_cd = build_coldata(picta_dna)
picta_rna_cd = build_coldata(picta_rna)
picta_dres = run_mf_deseq(picta_dna,
                          picta_dna_cd)
picta_rres = run_mf_deseq(picta_rna,
                          picta_rna_cd)
pd_res = format_res(picta_dres, 'picta','DNA')
pr_res = format_res(picta_rres, 'picta','RNA')
picta_fcs = rbind(pd_res,
                  pr_res)
picta_fcs %>% 
  mutate(sex_chrom = if_else(Chr==8,
                             'sex',
                             'auto')) %>% 
  ggplot(aes(x=sex_chrom, y=log2FoldChange, fill=source)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty=2) +
  # geom_hline(yintercept = c(-1,1)) +
  # coord_cartesian(y=c(-2,2)) +
  labs(y=bquote(log[2]~'M:F coverage'),
       subtitle = n) +
  coord_cartesian(y=c(-2, 2)) +
  theme(plot.subtitle = element_text(face = 'italic'),
        axis.title.x = element_blank())

#wrap them all
get_fc_res = function(spp){
  spp_dna = sub_counts(spp, 'DNA')
  spp_rna = sub_counts(spp, 'RNA')
  spp_dna_cd = build_coldata(spp_dna)
  spp_rna_cd = build_coldata(spp_rna)
  print('DNA codata:')
  print(spp_dna_cd)
  print('RNA coldata:')
  print(spp_rna_cd)
  spp_dres = run_mf_deseq(spp_dna,
                            spp_dna_cd)
  spp_rres = run_mf_deseq(spp_rna,
                            spp_rna_cd)
  d_res = format_res(spp_dres, spp,'DNA')
  r_res = format_res(spp_rres, spp,'RNA')
  spp_fcs = rbind(d_res,
                  r_res)
  return(spp_fcs)
}

plot_fc_boxes = function(fc_dat, spp){
  fc_dat %>% 
    mutate(sex_chrom = if_else(Chr==8,
                               'sex',
                               'auto')) %>% 
    ggplot(aes(x=sex_chrom, y=log2FoldChange, fill=source)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, lty=2) +
    labs(y=bquote(log[2]~'M:F coverage'),
         subtitle = spp) +
    coord_cartesian(y=c(-2, 2)) +
    theme(plot.subtitle = element_text(face = 'italic'),
          axis.title.x = element_blank())
}


#NOW RUN FOR RETICULATA AND WINGEI
picta_fcs
ret_fcs = get_fc_res('reticulata')
wing_fcs = get_fc_res('wingei')
pplt = plot_fc_boxes(picta_fcs, 'picta')
rplt = plot_fc_boxes(ret_fcs, 'reticulata')
wplt = plot_fc_boxes(wing_fcs, 'wingei')
plot_grid(rplt,wplt,pplt, nrow=3)



#LOAD THE PREVIOUS RNASEQ RESULTS
#reticulata
ll=load('rnaseq/results_files/PreticulataOnly_deseqResults.Rdata')
rres = format_res(res,
                  'reticulata',
                  'RNA')
ret_fcs2 = ret_fcs %>% 
  filter(source != 'RNA') %>% 
  rbind(rres) %>% 
  filter(!is.na(Chr))
#wingei
ll=load('rnaseq/results_files/PwingeiOnly_deseqResults.Rdata')
wres = format_res(res,
                  'wingei',
                  'RNA')


wing_fcs2 = wing_fcs %>% 
  filter(source != 'RNA') %>% 
  rbind(wres) %>% 
  filter(!is.na(Chr))
#picta
ll=load('rnaseq/results_files/MpictaOnly_deseqResults.Rdata')
pres = format_res(res,
                  'picta',
                  'RNA')
picta_fcs2 = picta_fcs %>% 
  filter(source != 'RNA') %>% 
  rbind(pres) %>% 
  filter(!is.na(Chr))

#PLOT WITH OLD RNA RESULTS
pplt = plot_fc_boxes(picta_fcs2, 'picta')
rplt = plot_fc_boxes(ret_fcs2, 'reticulata')
wplt = plot_fc_boxes(wing_fcs2, 'wingei')
plot_grid(rplt,wplt,pplt, nrow=3)
#These look right


# save the results for plotting selected figure ---------------------------

save(ret_fcs2, wing_fcs2, picta_fcs2, file='figure_plotting/fc_boxplot_objects.Rdata')

