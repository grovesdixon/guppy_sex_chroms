#!/usr/bin/env python

from sys import argv
import pandas as pd
import numpy as np
import subprocess
import os

#set up usage information
usage = "parse_vcf_depths.py infile.vcf outfile.tsv\n"
try:
    argv[1]
except IndexError:
    print('\nUsage:')
    exit(usage)

infile_name = argv[1]
outfile_name = argv[2]

#build temp file with header fixed for reading
print('Building pandas-readable version...')
temp_dir = 'tempdir_parse_vcf_depths'
if not os.path.isdir(temp_dir):
	os.mkdir(temp_dir)
temp_out = "{}/TEMP_{}_TEMP".format(temp_dir, infile_name)
with open(infile_name, 'r') as infile:
	with open(temp_out, 'w') as out:
		for line in infile:
			if line[0:2]=="#C":
				line = line.replace('#CHROM', 'CHROM')
			out.write(line)

#read in
print('Reading in...')
vdat = pd.read_csv(temp_out, sep='\t', comment='#')
cols = vdat.columns
geno_cols = cols[(list(cols).index('FORMAT')+1):]
gdat = vdat[geno_cols]

#get the AD index
print('parsing the allele depth index...')
ad_index = int(np.unique(vdat['FORMAT'].apply(lambda x: x.split(":").index('AD'))))

#pull out the depths for the reference and alternative alleles
print('parsing out the depths...')
def pull_ref_depth(x):
    return x.split(":")[ad_index].split(',')[0]
def pull_alt_depth(x):
    return x.split(":")[ad_index].split(',')[1]
ref_dat = gdat.applymap(pull_ref_depth)
alt_dat = gdat.applymap(pull_alt_depth)
ref_dat.columns = ['ref' + '_' + x for x in ref_dat.columns]
alt_dat.columns = ['alt' + '_' + x for x in alt_dat.columns]

#format output
print('writing results to {}...'.format(outfile_name))
pdat = vdat[['CHROM', 'POS', 'REF', 'ALT']]
out = pd.concat([pdat, ref_dat, alt_dat], axis=1)
out.to_csv(outfile_name, sep='\t', index=False)


