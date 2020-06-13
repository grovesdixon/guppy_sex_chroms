library(DESeq2)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

# EXPLORE -----------------------------------------------------------------

select_spp = 'reticulata'
select_spp = 'picta'
select_spp = 'wingei'

#LOAD
infile = paste('rnaseq/allele_specific_expression/', select_spp, '_xy_counts.Rdata', sep='')
ll=load(infile)
ll


head(summed_counts)

gene_means = apply(summed_counts, 1, mean)
keep = gene_means > 3
table(keep)


male_coldata = coldata[coldata$sex == 'male',]
male_counts = summed_counts[,rownames(male_coldata)]
gene_means = apply(male_counts, 1, mean)
keep = gene_means > 3
table(keep)
male_counts = male_counts[keep,]

#get variance stabilized counts for timepoint 2
dds<-DESeqDataSetFromMatrix(male_counts,
                            colData = male_coldata, 
                            design = formula(~ chrom + individual)
)

INDEPDENDENT_FILTERING = FALSE
dds <- DESeq(dds)
resultsNames(dds)
res = results(dds, contrast = c('chrom', 'Y', 'X'), independentFiltering=INDEPDENDENT_FILTERING)
rld = rlog(dds)

#volcano for tissue
TITLE='Y vs X'
vol = data.frame(res) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  ggplot(aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1) +
  scale_color_manual(values=c('black', 'red')) + 
  labs(subtitle='Y vs X',
       y=bquote(-log[10]~pvalue),
       x=bquote(log[2]~fold~difference))  +
  theme(legend.position='none')
vol


#LOAD POSITION DATA
gene_coords = read_tsv('metadata/geneCoords.tsv',
                       col_names = c('chr', 'start', 'end', 'gene'))

res_df = data.frame(res) %>% 
  rownames_to_column('gene') %>% 
  mutate(gene = sub('gene:', '', gene)) %>% 
  as_tibble()

print('all genes have coords?')
print(sum(res_df$gene %in% gene_coords$gene) == nrow(res_df))

res_df = res_df %>% 
  left_join(gene_coords, by = 'gene')


res_df %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  ggplot(aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1) +
  scale_color_manual(values=c('black', 'red')) + 
  labs(subtitle='Y vs X',
       y=bquote(-log[10]~pvalue),
       x=bquote(log[2]~fold~difference))  +
  theme(legend.position='none') +
  facet_wrap(~chr)


#DENSITY PLOT OF Y VS X
#for all chroms
res_df %>% 
  ggplot(aes_string(x='log2FoldChange')) +
  geom_histogram() +
  geom_vline(xintercept = 0, lty=2) +
  facet_wrap(~chr)

#by chrom type
res_df %>% 
  mutate(sex_chrom = if_else(chr==8,
                             'sex',
                             'autosome')) %>% 
  ggplot(aes_string(x='log2FoldChange', fill='sex_chrom')) +
  geom_density(alpha=0.75) +
  geom_vline(xintercept = 0, lty=2) +
  scale_fill_manual(values=c('grey50', 'goldenrod')) +
  labs(x=bquote('log'[2]~'fold difference'),
       fill='')


#CHECK THE Y-LINKED CALLING SUMMARY

chr_lengths = read_tsv('metadata/chrLengths.tsv',
                       col_names=c('chr', 'length'))

ysum_infile = paste('rnaseq/callXY/', select_spp, '_ylike_summary.tsv', sep='')
ysum = read_tsv(ysum_infile)
ysum %>% 
  filter(grepl('male', source)) %>% 
  mutate(chr=as.character(chr)) %>% 
  left_join(chr_lengths, by = 'chr') %>% 
  mutate(ylike_density = N_ylike / length,
         chr = factor(chr, levels = as.character(1:24))) %>% 
  ggplot(aes(x=chr, y = ylike_density)) +
  geom_bar(stat='identity') +
  labs(y='density of Y-like loci',
       x='chromosome')

#CHECK WHERE THE Y-LIKE ARE
ylike_in = paste('rnaseq/allele_specific_expression/all_SNPs_', select_spp, '_y_X_maculatus5.0.txt', sep='')
ylike = read_tsv(ylike_in,
                 col_names = c('id', 'chr', 'pos', 'xx', 'X.Y')) %>% 
  separate(X.Y,
           into = c('X', 'Y'))

#accross all chroms
ylike %>% 
  ggplot(aes(x=pos/1e6)) +
  geom_histogram() +
  labs(y='Y-like locus count',
       x='position (Mb)') +
  facet_wrap(~chr)

#for just the sex chromosome
ylike %>% 
  filter(chr==8) %>% 
  mutate(revMb = (max(pos)-pos) / 1e6) %>% 
  ggplot(aes(x=revMb)) +
  geom_histogram() +
  labs(y='Y-like locus count',
       x='position (Mb)') +
  facet_wrap(~chr)




