#plot_pcas.R

rm(list=ls())
source('functions.R')


#GLOBAL VARS
NTOP = 10000
FIX_COORDS = FALSE


# LOAD THE RNA RLD DATA ---------------------------------------------------

#PLOT FOR ALL SAMPLES AT ONCE
library(DESeq2)
ll=load('rnaseq/results_files/all_rld.Rdata')
ll
pca_df = build_pca(rld.df, coldata, ntop = NTOP, pcs = 2)
addx=3
addy=2
percentVar = attr(pca_df, "percentVar")[c(1,2)]
xlab = paste0(paste0(paste0("PC", 1), ": "), round(percentVar[1] * 100), "%")
ylab = paste0(paste0(paste0("PC", 2), ": "), round(percentVar[2] * 100), "%")
pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=Organism, shape=sex)) +
  geom_point(size=5) +
  lims(x=c(min(pca_df$PC1)-addx, max(pca_df$PC1+addx)),
       y=c(min(pca_df$PC2)-addy, max(pca_df$PC2+addy))) +
  labs(x=xlab, y=ylab)


#LOAD DATA FOR EACH SPECIES INDIVIDUALLY
dat_names = c("PreticulataOnly",
              "PwingeiOnly",
              "MpictaOnly",
              "PlatipinnaOnly",
              "GholbrookiOnly") #reorder them to match figures
path_list = paste('rnaseq/results_files/', dat_names, '_rld.Rdata', sep='')
names(path_list) = dat_names

#funciton to read in rld data
load_rld = function(x){
  print(paste('loading ', x, '...', sep=''))
  ll=load(x)
  rld.df = data.frame(assay(rld))
  print('Names align?')
  print(sum(colnames(rld.df) == rownames(coldata))==ncol(rld.df))
  return(list('rld'=rld.df,
              'coldata'=coldata))
}

#get the rld data for each project
rld_list = map(path_list, load_rld)

#funnction to plot the PCA
rna_pca_list = list()
for (n in dat_names){
  x=rld_list[[n]]
  rld = x[['rld']]
  coldata = x[['coldata']]
  pca_df = build_pca(rld, coldata, ntop = NTOP, pcs = 2)
  plt = plot_rld_pca(pca_df,
                     group_col = 'sex',
                     fix_coords = FIX_COORDS) + 
    theme(legend.position = 'right') +
    labs(subtitle = n)
  rna_pca_list[[n]] = plt
}

plot_grid(plotlist = rna_pca_list)
#re-assuringly (at least for reticulata, wingei, and picta), males and females split on PC1



#check the weird picta sample
picta = rld_list[["MpictaOnly"]]
coldata = data.frame(picta['coldata'])
coldata$sample = rownames(coldata)
pca_df = build_pca(picta[['rld']], coldata, ntop = NTOP, pcs = 2)
pca_df %>% 
  filter(PC2 == min(PC2))




