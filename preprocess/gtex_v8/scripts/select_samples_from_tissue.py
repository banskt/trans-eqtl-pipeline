#!/usr/bin/env python3
# Author: Saikat Banerjee

import pandas as pd
import argparse
import gzip
import re

def parse_args():

    parser = argparse.ArgumentParser(description='Extract samples of a particular tissue from GTEx gene expression file')

    parser.add_argument('--rpkm',
                        type=str,
                        dest='rpkmpath',
                        metavar='FILE',
                        help='input GCT file containing all samples (RPKMs)')

    parser.add_argument('--counts',
                        type=str,
                        dest='countspath',
                        metavar='FILE',
                        help='input GCT file containing all samples (Read Counts)')

    parser.add_argument('--output',
                        type=str,
                        dest='outfilepath',
                        metavar='FILE',
                        help='output GCT file containing selected samples')

    parser.add_argument('--tissue',
                        type=str,
                        dest='tissue',
                        metavar='STR',
                        help='exact name of the tissue')

    parser.add_argument('--pheno',
                        type=str,
                        dest='phenofilepath',
                        metavar='STR',
                        help='phenotype file of the GTEx consortium')

    opts = parser.parse_args()
    return opts


def get_samples(filepath, tissue):
    """
    Extract sample IDs from GTEx phenotype description
    """
    with gzip.open(filepath, 'r') as mfile:
        df = pd.read_csv(filepath, sep='\t', skiprows=10, index_col=0)
        df = df.loc[(df['SMTORMVE'] != "FLAGGED") & (df['SMGEBTCHT'] == "TruSeq.v1") & (df['SMTSD'] == tissue)]
    samples = df['SAMPID'].tolist()
    return samples


def read_gct(gct_file, sample_ids):
    """
    Load GCT as DataFrame. Keep the description, and only the samples of sample_id.
    """
    print("Reading ", gct_file)
    if re.search("gz", gct_file):
        ins = gzip.open(gct_file)
        skip = 0
        for line in ins:
            if re.search("^name", line.decode(), re.IGNORECASE): break
            else: skip += 1
        ins.close()
    else:
        ins = open(gct_file)
        skip = 0
        for line in ins:
            if re.search("^name", line, re.IGNORECASE): break
            else: skip += 1
        ins.close()
    
    df = pd.read_csv(gct_file, sep='\t', skiprows=skip, index_col=0)
    if "Description" in df.columns:
        df.drop('Description', axis=1, inplace=True)
    df.index.name = 'gene_id'
    
    # select the samples
    df = df[[i for i in df.columns if i in sample_ids]]
    
    # add in the descriptions and rearrange the columns
    # df['Description'] = description_df
    # cols = df.columns.tolist()  
    # newcols = cols[-1:] + cols[:-1]
    return df #[newcols]

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

    sample_ids = get_samples(opts.phenofilepath, opts.tissue)
    print ('Number of unflagged, RNA-Seq samples from %s: %i\n' % (opts.tissue, len(sample_ids)))
    df = read_gct(opts.rpkmpath, sample_ids)
    print ('Number of samples read from GCT file (RPKMs): %i\n' % (df.shape[1] - 1))
    write_gct(df, opts.outfilepath+"_tpm.gct")
    print ('New GCT file written with %i samples and %i genes.\n' % (df.shape[1], df.shape[0] - 1))

    df = read_gct(opts.countspath, sample_ids)
    print ('Number of samples read from GCT file (Read Counts): %i\n' % (df.shape[1] - 1))
    write_gct(df, opts.outfilepath+"_counts.gct")
    print ('New GCT file written with %i samples and %i genes.\n' % (df.shape[1], df.shape[0] - 1))

