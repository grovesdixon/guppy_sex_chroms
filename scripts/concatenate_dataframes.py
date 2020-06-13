#!/usr/bin/env python
##concatenate_dataframes.py
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
    Take a regular expression for a set of files and concatenate them into single dataframe

    '''

    parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
    parser.add_argument('-i', required = True, nargs = '+', dest = 'infile_names', help = 'Regular expression for multiple files. They should have the same header')
    parser.add_argument('-o', required = True, dest = 'outfile_name', help = 'The the input file')

    #--- PARSE ARGUMENTS ---#
    args = parser.parse_args()
    infile_list = args.infile_names
    outfile_name = args.outfile_name

    print('\nConcatentating datframes from the following files:')
    for i in infile_list:
        print(i)

    first = infile_list[0]
    dat = pd.read_csv(first, sep='\t')
    rest = infile_list[1:]
    for f in rest:
        new_dat = pd.read_csv(f, sep='\t')
        if np.sum(new_dat.columns==dat.columns)==len(dat.columns):
            dat = pd.concat([dat, new_dat], axis=0)
        else:
            print('\nERROR. Columns to do match between tables.')
            print('check table {}'.format(f))
            exit()

    print('\nGood, all files had same columns.')
    print('writing results to {}...'.format(outfile_name))
    dat.to_csv(outfile_name, sep='\t', index=False)
    print('done')













