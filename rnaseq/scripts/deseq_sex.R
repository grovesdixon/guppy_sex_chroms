#deseq_sex.R 

#SET UP THE DATA TO RUN DESEQ
library(DESeq2)
library(cowplot)
library(tidyverse)
rm(list=ls())

source('rnaseq/scripts/rnaseq_functions.R')

# LOAD THE SUBSET DO ANALYZE  -----------------------------------------------------------
#(see rnaseq/sex_compare_groups for these)
source('rnaseq/sex_compare_groups/Gholbrooki.R')
source('rnaseq/sex_compare_groups/Mpicta.R') 
source('rnaseq/sex_compare_groups/Platipinna.R') 
source('rnaseq/sex_compare_groups/Preticulata.R') 
source('rnaseq/sex_compare_groups/Pwingei.R') 

# RUN DESEQ ---------------------------------------------------------------

coldata$Organism = factor(coldata$Organism)
coldata$sex = factor(coldata$sex)

#get variance stabilized counts for timepoint 2
dds<-DESeqDataSetFromMatrix(counts,
	colData = coldata, 
	design = formula(~ sex)
	)

INDEPDENDENT_FILTERING = FALSE
dds <- DESeq(dds)
resultsNames(dds)
res = results(dds, contrast = c('sex', 'male', 'female'), independentFiltering=INDEPDENDENT_FILTERING)
rld = rlog(dds)

#volcano for tissue
vol = data.frame(res) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  ggplot(aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
	geom_point(alpha=0.1) +
	scale_color_manual(values=c('black', 'red')) + 
	labs(title = TITLE,
	     subtitle='M vs F',
	     y=bquote(-log[10]~pvalue),
	     x=bquote(log[2]~fold~difference))  +
  theme(legend.position='none')


#save results
outname = paste('rnaseq/results_files/', PREFIX, '_deseqResults.Rdata', sep='')
save(res,coldata, file=outname)
rldOut = paste('rnaseq/results_files/', PREFIX, '_rld.Rdata', sep='')
save(rld,coldata, file=rldOut)



