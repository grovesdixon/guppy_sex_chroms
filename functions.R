#functions.R


# LIBS --------------------------------------------------------------------

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())



# GLOBAL VARS -------------------------------------------------------------

species_names = c('reticulata',
                  'wingei',
                  'picta',
                  'latipinna',
                  'Gambusia')

sex_col_list = list('sex' = 'goldenrod',
                    'autosome' = 'grey50')


# FUNCTIONS ---------------------------------------------------------------

#function to convert the hellerii Genebank chromosome tags to numbers
swap_chrs = function(df){
  ll=load('metadata/hellerii_chr_tags.Rdata')
  df %>% 
    mutate(chr=hellerii_chr_tags[df$chr])
}


rem_legend = function(x){
  return(x + theme(legend.position = 'none'))
}


get_95_conf_intervals = function(vector){
  q = quantile(vector, probs=seq(0,1,0.05), na.rm=TRUE)
  lower = q['5%']
  upper = q['95%']
  return(c(lower, upper))
}


get_95_conf_intervals_two_tailed = function(vector){
  q = quantile(vector, probs=seq(0,1,0.005), na.rm=TRUE)
  lower = q['2.5%']
  upper = q['97.5%']
  return(c(lower, upper))
}


#ass lower and upper columns to a dataframe based on a factor (here usually species)
#idea is to get an upper and lower bound for a given factor level in a dataframe,
#eg the 95% confidence intervals for FST for each species in a df
add_intervals_by_factor_ALL = function(dat, factor_col, value_col){
  factors = dat %>% 
    pull(factor_col)
  lowers = c()
  uppers = c()
  for (f in unique(factors)){
    v = dat[factors==f,value_col]
    q = get_95_conf_intervals(v)
    lowers[f] = q[1]
    uppers[f] = q[2]
  }
  dat %>% 
    mutate(lower = lowers[factors],
           upper = uppers[factors])
}

#sampe as above but uses all but the sex chromosome for confidence intervals
add_intervals_by_factor = function(dat, factor_col, value_col, sexchr=8){
  factors = dat %>% 
    pull(factor_col)
  lowers = c()
  uppers = c()
  for (f in unique(factors)){
    v = dat[factors==f & dat$chr!=sexchr, value_col]
    q = get_95_conf_intervals(v)
    lowers[f] = q[1]
    uppers[f] = q[2]
  }
  dat %>% 
    mutate(lower = lowers[factors],
           upper = uppers[factors])
}


#sampe as above but uses all but the sex chromosome for confidence intervals
add_intervals_by_factor_two_tailed = function(dat, factor_col, value_col, sexchr=8){
  factors = dat %>% 
    pull(factor_col)
  lowers = c()
  uppers = c()
  for (f in unique(factors)){
    v = dat[factors==f & dat$chr!=sexchr, value_col]
    q = get_95_conf_intervals_two_tailed(v)
    lowers[f] = q[1]
    uppers[f] = q[2]
  }
  dat %>% 
    mutate(lower = lowers[factors],
           upper = uppers[factors])
}


#function to plot the sex chromosome scatterplot with 95 confidence intervals
plot_sexchrom_scatter = function(df0, xcol, ycol, factor_col='species', sexchrom=8, ylim=FALSE,alpha=1){
  df = add_intervals_by_factor(df0, factor_col, ycol)
  y = df %>% 
    pull(ycol)
  df$upper_outlier = y > df$upper
  plt = df %>% 
    filter(chr==sexchrom) %>% 
    ggplot() +
    # geom_hline(yintercept = 0, lty=2) +
    geom_ribbon(aes_string(x=xcol, ymin='lower', ymax='upper'),
                fill='grey',
                color='grey') +
    geom_point(aes_string(x=xcol, y=ycol, color='upper_outlier'), alpha=alpha) +
    geom_smooth(aes_string(x=xcol, y=ycol), lwd=0.75, se=FALSE) +
    scale_color_manual(values=c('black', 'black')) +
    theme(legend.position = 'none')
  if (length(ylim) > 1){
    plt=plt+lims(y=ylim)
  }
  return(plt)
}

#same as above but for two-tailed
plot_sexchrom_scatter_two_tailed = function(df0, xcol, ycol, factor_col='species', sexchrom=8, ylim=FALSE,alpha=1){
  df = add_intervals_by_factor_two_tailed(df0, factor_col, ycol)
  y = df %>% 
    pull(ycol)
  df$outlier = y > df$upper | y < df$lower
  plt = df %>% 
    filter(chr==sexchrom) %>% 
    ggplot() +
    geom_ribbon(aes_string(x=xcol, ymin='lower', ymax='upper'), fill='grey', color='grey') +
    # geom_hline(yintercept = 0, lty=2) +
    geom_point(aes_string(x=xcol, y=ycol, color='outlier'), alpha=alpha) +
    geom_smooth(aes_string(x=xcol, y=ycol), lwd=0.75, se=FALSE) +
    scale_color_manual(values=c('black', 'black')) +
    theme(legend.position = 'none')
  if (length(ylim) > 1){
    plt=plt+lims(y=ylim)
  }
  return(plt)
}


#function to plot pca from normalized expression counts
build_pca = function(df,
                     coldata,
                     ntop = 25000,
                     pcs = 6){
  #get row varainces and select top
  df = as.matrix(df)
  rv <- rowVars(df)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  #build pca
  pca <- prcomp(t(df[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  d = cbind(data.frame(pca$x[,1:pcs]),
            coldata)
  attr(d, "percentVar") <- percentVar[1:2]
  return(d)
}



#function to plot PCA scatterplot from output from build_rld_pca
plot_rld_pca = function(df,
                        group_col = 'treatment',
                        pc1 = 1,
                        pc2 = 2,
                        subtitle = "",
                        size = 4,
                        legend_title=NULL,
                        x_invert=1,
                        legend_position = 'none',
                        fix_coords = TRUE){
  #select PCs to plot and nvert X if needed
  plt_df = tibble(x = df[,paste('PC', pc1, sep = '')]*x_invert,
                  y = df[,paste('PC', pc2, sep = '')],
                  col = factor(df[,group_col]))
  #pull out the percent variances
  percentVar = attr(df, "percentVar")[c(pc1,pc2)]
  #build axis labs
  xlab = paste0(paste0(paste0("PC", pc1), ": "), 
                round(percentVar[pc1] * 100), "%")
  ylab = paste0(paste0(paste0("PC", pc2), ": "), round(percentVar[pc2] * 100), "%")
  g = plt_df %>% 
    ggplot(aes(x=x,
               y=y,
               color=col)) + 
    geom_point(size = size) +
    labs(x = xlab,
         y = ylab,
         subtitle=subtitle,
         color=legend_title) +
    theme(legend.position=legend_position,
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank())
  if (fix_coords){
    g = g + coord_fixed() 
  }
  return(g)
}
