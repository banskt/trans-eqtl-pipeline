#!/usr/bin/env python3
# Author: Saikat Banerjee

import pandas as pd
import argparse
import gzip
import os
import numpy as np

def parse_args():

    parser = argparse.ArgumentParser(description='Extract samples of a particular tissue from GTEx gene expression file')

    parser.add_argument('--tpm',
                        type=str,
                        dest='tpmfilepath',
                        metavar='FILE',
                        help='input GCT file containing TPM data')

    parser.add_argument('--counts',
                        type=str,
                        dest='countsfilepath',
                        metavar='FILE',
                        help='input GCT file containing counts data')

    parser.add_argument('--outdir',
                        type=str,
                        dest='outdir',
                        metavar='STR',
                        help='path of the output directory')

    opts = parser.parse_args()
    return opts


def read_gct(gct_file):
    """
    Load GCT as DataFrame.
    """
    df = pd.read_csv(gct_file, sep='\t', skiprows=2, index_col=0)
    description_df = df['Description']
    df.index.name = 'gene_id'
    
    # add in the descriptions and rearrange the columns
    #df['Description'] = description_df
    #cols = df.columns.tolist()  
    #newcols = cols[-1:] + cols[:-1]

    # remove descriptions
    cols = df.columns.tolist()
    newcols = cols[1:]

    return df[newcols]

def write_gct(df, filepath):
    """
    Write dataframe as a GCT file
    """
    with open(filepath, 'w') as mfile:
        mfile.write("#1.2\n")
        mfile.write('%i\t%i\n' % (df.shape[0], df.shape[1] - 1))
        # mfile.write(str(df.shape[0])+'\t'+str(df.shape[1] - 1)+'\n')
        df.to_csv(mfile, sep='\t', index=True, header=True)
    

if __name__=='__main__':

    opts = parse_args()
    min_samples = 0
    count_threshold = 0

    # Discard genes with all-zero counts
    print ("Reading GCT files ...")
    tpm_df = read_gct(opts.tpmfilepath)
    counts_df = read_gct(opts.countsfilepath)
    print ("Completed reading. Writing new files ...")
    mask = (np.sum(counts_df > count_threshold, axis=1) > min_samples).values

    counts_outfile = os.path.join(opts.outdir, "counts_nozero.gct")
    tpm_outfile = os.path.join(opts.outdir, "tpm_nozero.gct")

    write_gct(counts_df[mask], counts_outfile)
    write_gct(tpm_df[mask], tpm_outfile)
    print ('New GCT file written with %i samples and %i genes.\n' % (counts_df[mask].shape[1], counts_df[mask].shape[0] - 1))
