#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 16:21:37 2020

@author: grovesdixon
"""


import pandas as pd
import numpy as np
import os
os.chdir('/Users/grovesdixon/gitreps/guppy_sex_chroms/')

infile_name = '/Users/grovesdixon/Desktop/parse_depth_vcf/chr8_all_wingei_RNA_depths.tsv'
gene_coords_in = 'metadata/geneCoords.tsv'
outfile_name = '/Users/grovesdixon/Desktop/parse_depth_vcf/chr8_all_wingei_RNA_depthsByGene.tsv'


#read in the depth data (use parse_vcf_depths.py to build this)
ddat = pd.read_csv(infile_name, sep='\t')
ddat['CHROM'] = ddat['CHROM'].astype('str')
pdat = ddat[['CHROM', 'POS', 'REF', 'ALT']]
gdat = ddat[ddat.columns[5:]]


#read in the gene coordinates
cdat = pd.read_csv(gene_coords_in, sep='\t',
                   names = ['CHROM', 'start', 'end', 'name'])

#subset by chromosome (in case the depth file isn't already)
uchrs = np.unique(ddat['CHROM'])
print('chromosomes found in the input file:')
print(uchrs)

#get means for each gene
cols = list(ddat.columns[4:])
cols.append('gene')
mn_dat = pd.DataFrame(columns=cols)

for CHR in uchrs:
    print('assigning genes from chr {}'.format(CHR))
    chr_cdat = cdat.loc[cdat['CHROM']==CHR]
    chr_ddat = ddat.loc[ddat['CHROM']==CHR]
    for i in chr_cdat.index.values:
        row = chr_cdat.loc[i]
        sub = chr_ddat.loc[(chr_ddat['POS']>=row['start']) & (chr_ddat['POS'] <= row['end'])]
        sub_depth = sub[sub.columns[4:]]
        mn_depth = sub_depth.mean(axis=0)
        mn_depth = pd.DataFrame([mn_depth])
        mn_depth[mn_depth.isnull()] = np.NaN
        mn_depth['gene'] = row['name']
        mn_dat = mn_dat.append(mn_depth)
        
mn_dat.head()
mn_dat.to_csv(outfile_name, sep='\t', na_rep='NA', index=False)

