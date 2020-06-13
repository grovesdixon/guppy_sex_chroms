#MvF_along_chroms.R
source('functions.R')

# LOAD GENE COORDINATES -------------------------------------------------------------
ldat = read_tsv('metadata/chrLengths.tsv',
                col_names = c('chr', 'length')) %>% 
  filter(chr %in% 1:24)
lengths = ldat$length
names(lengths) = ldat$chr
gdat = read_tsv('metadata/geneCoords.tsv',
                col_names=c('chr', 'start', 'end', 'gene')) %>% 
  filter(chr %in% 1:24)
gdat$chr=as.character(gdat$chr)

uchrs = na.omit(unique(gdat$chr))
gdat$length = 0
for (ichr in uchrs){
  l=lengths[ichr]
  gdat = gdat %>% 
    mutate(length = if_else(chr==ichr,
                            l,
                            length))
}
gdat$scaled = gdat$start / gdat$length


# LOAD DESEQ RESULTS ------------------------------------------------------

prefixes = c('GholbrookiOnly',
             'MpictaOnly',
             'PlatipinnaOnly',
             'PreticulataOnly',
             'PwingeiOnly')

rldList = list()
#populate rldList
for (p in prefixes){
  ll=load(paste('rnaseq/results_files/',p,'_rld.Rdata',sep=''))
  
  #load, format, and add rld data
  rld.df=assay(rld) %>% 
    data.frame() %>% 
    mutate(gene=rownames(assay(rld))) %>% 
    left_join(gdat, by = 'gene') %>% 
    as_tibble()
  rldList[[p]] = rld.df
}

more.prefixes = c('GholbrookiOnly',
                  'MpictaOnly',
                  'PlatipinnaOnly',
                  'PreticulataOnly',
                  'PwingeiOnly')

resList = list()
#poulate resList
for (p in more.prefixes){
  ll=load(paste('rnaseq/results_files/',p,'_deseqResults.Rdata',sep=''))
  
  #load, and format res data
  res.df = res %>% 
    data.frame() %>% 
    mutate(gene=rownames(res)) %>% 
    left_join(gdat, by='gene') %>% 
    filter(!is.na(chr)) %>% 
    as_tibble()
  chrNums = as.numeric(res.df$chr)
  lastChr = max(chrNums, na.rm=TRUE)
  res.df$chr=factor(res.df$chr, levels=1:lastChr)
  
  #designate percentiles for log2fold
  q=quantile(res.df$log2FoldChange, probs=seq(0,1,0.025), na.rm=TRUE)
  low = q['2.5%']
  high = q['97.5%']
  res.df$extremeL2 = res.df$log2FoldChange < low | res.df$log2FoldChange > high
  #designate for abs log2fold
  res.df$abslog2 = abs(res.df$log2FoldChange)
  absq = quantile(res.df$abslog2, probs=seq(0,1,0.05), na.rm=TRUE)
  abshigh = absq['95%']
  res.df$extremeAbsL2 = res.df$log2FoldChange < low | res.df$log2FoldChange > high
  
  #add to list
  resList[[p]] = res.df
}


# SAVE --------------------------------------------------------------------

tags = names(resList)
spp_names = c('Gambusia',
              'picta',
              'latipinna',
              'reticulata',
              'wingei')
mod_res = list()
mod_rld = list()
for (i in 1:length(tags)){
  tag = tags[i]
  sn = spp_names[i]
  print('------')
  print(tag)
  print(sn)
  rld = rldList[[i]]
  res = resList[[i]]
  rld$species = sn
  res$species = sn
  mod_rld[[sn]] = rld
  mod_res[[sn]] = res
}

ge_res = purrr::reduce(mod_res, rbind)
ge_rld = purrr::reduce(mod_rld, full_join, by=c('gene', 'chr', 'start', 'end', 'length', 'scaled'))
save(ge_res, file='figure_plotting/ge_res.Rdata')
save(ge_rld, file='figure_plotting/ge_rld.Rdata')

# PLOT P. PICTA SEX DIFFERENCES ---------------------------------------


p='MpictaOnly'
rld.df = rldList[[p]]
res.df = resList[[p]]

#boxplot log2foldchanges
res.df %>% 
  mutate(sex = if_else(chr=='8',
                       'sex',
                       'autosome')) %>% 
  ggplot(aes(x=chr, y=log2FoldChange, fill=sex)) +
  geom_hline(yintercept = 0, lwd=0.5, lty=3, color='black') +
  geom_boxplot() +
  labs(x='chromosome',
       y=bquote("M:F"~log[2]~"fold change"),
       title = 'Differential expression between sexes',
       subtitle='P. picta') +
  lims(y=c(-1,1)) +
  theme(plot.subtitle = element_text(face ='italic'))


#boxplot abs log2foldchanges
res.df %>% 
  ggplot(aes(x=chr, y=abslog2)) +
  geom_boxplot()


#plot along sex chromosome
res.df %>% 
  filter(!is.na(log2FoldChange),
         chr==8) %>% 
  ggplot(aes(x=start,y=log2FoldChange)) +
  geom_hline(yintercept = 0, lty=2) +
  geom_smooth(se=FALSE, method = 'loess', span=0.1)+
  geom_point(aes(color=extremeL2)) +
  labs(x='Chr 8 position',
       y=bquote("M:F"~log[2]~"fold change"),
       title = 'Differential expression between sexes',
       subtitle='P. picta',
       color = '95 percentile') +
  theme(plot.subtitle = element_text(face ='italic'))
  


# PLOT PICTA VS Gholbrooki ------------------------------------------------

#function to load the rld data and get the gene means by Organism
load_rld = function(rldPath, outSpp, compareSpp){
  ll=load(rldPath)
  rld.df = data.frame(assay(rld))
  out = rownames(coldata)[coldata$Organism==outSpp]
  compare = rownames(coldata)[coldata$Organism==compareSpp]
  print('Outgroup samples:')
  print(out)
  print('Compare group samples:')
  print(compare)
  
  #convert to long format
  rld.long = rld.df %>% 
    mutate(gene=rownames(rld.df)) %>% 
    pivot_longer(cols=1:ncol(rld.df),
                 names_to='run',
                 values_to = 'lvl') %>% 
    mutate(Organism = if_else(run %in% out,
                              outSpp,
                              compareSpp))
  
  #get means
  mns = rld.long %>% 
    group_by(gene, Organism) %>% 
    summarize(mn = mean(lvl, na.rm=TRUE)) %>% 
    left_join(gdat, by = 'gene') %>% 
    filter(!is.na(chr)) %>% 
    pivot_wider(id_cols=c(gene, Organism, chr, start, end, length, scaled),
                names_from = Organism,
                values_from = mn)
}

#get female mean rlds
f.lvls = load_rld(rldPath = 'rnaseq/results_files/Ppicta_v_Gholbrooki_females_rld.Rdata',
                 outSpp = 'Gambusia_holbrooki',
                 compareSpp = 'Micropoecilia_picta')

#get male mean rlds
m.lvls = load_rld(rldPath = 'rnaseq/results_files/Ppicta_v_Gholbrooki_males_rld.Rdata',
                 outSpp = 'Gambusia_holbrooki',
                 compareSpp = 'Micropoecilia_picta')


#function to plot the level comparison
plot_lvl_compare = function(df.lvl, sexChrom, outSpp, compareSpp, xlab, ylab, subtitle){
  df.lvl$sex = if_else(df.lvl$chr == sexChrom,
                       'sex',
                       'autosome')
  sex.df = df.lvl %>% 
    filter(sex=='sex')
  plt = df.lvl %>% 
    ggplot(aes_string(x=outSpp, y=compareSpp, color='sex')) +
    geom_point() +
    scale_color_manual(values=c('black', 'red')) +
    labs(x=xlab,
         y=ylab,
         subtitle=subtitle) +
    theme(axis.title = element_text(face='italic'))
  plt = plt +
    geom_point(data=sex.df, aes_string(x=outSpp, y=compareSpp), color='red') +
    geom_abline(slope=1, intercept=0, color='grey', lty=1, lwd=1)
  return(plt)
}

fscatter = plot_lvl_compare(f.lvls,
                 sexChrom='8',
                 outSpp = 'Gambusia_holbrooki',
                 compareSpp = 'Micropoecilia_picta',
                 xlab = 'G. holbrooki',
                 ylab = 'P. picta',
                 subtitle='female')

mscatter = plot_lvl_compare(m.lvls,
                 sexChrom='8',
                 outSpp = 'Gambusia_holbrooki',
                 compareSpp = 'Micropoecilia_picta',
                 xlab = 'G. holbrooki',
                 ylab = 'P. picta',
                 subtitle='male')


plot_grid(fscatter + theme(legend.position = 'none'),
          mscatter,
          rel_widths=c(0.8, 1))



#load male data
# m.rld.out = get_mn_lvl('rnaseq/results_files/Gambusia_holbrooki_male_rld.Rdata')
# m.rld.c = get_mn_lvl('rnaseq/results_files/Micropoecilia_picta_male_rld.Rdata')
ll = load('rnaseq/results_files/Ppicta_v_Gholbrooki_males_deseqResults.Rdata')
mres = data.frame(res) %>% 
  mutate(gene=rownames(res),
         sex='male') %>% 
  left_join(gdat, by = 'gene')


#load female data
# f.rld.out = get_mn_lvl('rnaseq/results_files/Gambusia_holbrooki_female_rld.Rdata')
# f.rld.c = get_mn_lvl('rnaseq/results_files/Micropoecilia_picta_female_rld.Rdata')
ll = load('rnaseq/results_files/Ppicta_v_Gholbrooki_females_deseqResults.Rdata')
fres = data.frame(res) %>% 
  mutate(gene=rownames(res),
         sex='female') %>% 
  left_join(gdat, by = 'gene')


#combine
res = rbind(mres, fres) %>% 
  filter(!is.na(chr))
res$sexChrom = if_else(res$chr=='8',
                       'sex',
                       'autosome')

res %>% 
  ggplot(aes(x=sexChrom, fill=sex, y=log2FoldChange)) +
  geom_hline(yintercept = 0, lty=2) +
  geom_boxplot() +
  lims(y=c(-3,3)) +
  labs(y=bquote(log[2]~'P.picta : G.holbrooki'),
       x='Chromosome',
       title='Comparison to ancestral state')
