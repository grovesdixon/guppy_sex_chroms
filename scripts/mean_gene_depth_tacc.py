#!/usr/bin/env python
##mean_gene_depth_tacc.py
##Groves Dixon

#import modules
import argparse
import pandas as pd
import numpy as np


##################################
############## MAIN ##############
##################################

if __name__ == '__main__':


#START RUN TIME CLOCK

    ##SET UP ARGUMENT PARSING
    Description = '''
    Description:
    Take a set of allele depths from a vcf (parsed out with parse_vcf_depths.py)
    and a set of gene coordinates and get the mean depth for each gene.

    '''

    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-i', required = True, dest = 'infile_name', help = 'The the input file (build this with parse_vcf_depths.py)')
    parser.add_argument('-g', required = True, dest = 'gene_coords', help = 'Table of gene coordiantes (bed file with chr\tstart\tend\tgeneName)')
    parser.add_argument('-o', required = True, dest = 'outfile_name', help = 'The the input file')

    #--- PARSE ARGUMENTS ---#
    args = parser.parse_args()
    infile_name = args.infile_name
    gene_coords_in = args.gene_coords
    outfile_name = args.outfile_name



    #---- RUN FUNCTIONS ----#
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











