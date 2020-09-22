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


# isolate positions and counts -------------------------------------------------------

pos = full_dat %>%
  dplyr::select(Geneid:Length)


counts = full_dat %>%
  dplyr::select(7:ncol(full_dat))

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
picta_dna_cd = build_coldata(picta_dna)
picta_dres = run_mf_deseq(picta_dna,
                          picta_dna_cd)
pd_res = format_res(picta_dres, 'picta','DNA')
picta_fcs = pd_res


#wrap them all
get_fc_res = function(spp){
  spp_dna = sub_counts(spp, 'DNA')
  # spp_rna = sub_counts(spp, 'RNA')
  spp_dna_cd = build_coldata(spp_dna)
  # spp_rna_cd = build_coldata(spp_rna)
  print('DNA codata:')
  print(spp_dna_cd)
  spp_dres = run_mf_deseq(spp_dna,
                            spp_dna_cd)
  d_res = format_res(spp_dres, spp,'DNA')
  return(d_res)
}

# #NOW RUN FOR RETICULATA AND WINGEI
picta_fcs
ret_fcs = get_fc_res('reticulata')
wing_fcs = get_fc_res('wingei')


#LOAD THE PREVIOUS RNASEQ RESULTS USED FOR SCATTERPLOTS
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

#PLOT RESULTS

#function to plot the boxplots
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

#build plots
pplt = plot_fc_boxes(picta_fcs2, 'picta')
rplt = plot_fc_boxes(ret_fcs2, 'reticulata')
wplt = plot_fc_boxes(wing_fcs2, 'wingei')
plot_grid(rplt,wplt,pplt, nrow=3)



# save the results for plotting selected figure ---------------------------

save(ret_fcs2, wing_fcs2, picta_fcs2, file='figure_plotting/fc_boxplot_objects.Rdata')

