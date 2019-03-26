import sys, os
import readgtf
from collections import defaultdict
import pandas as pd
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(description='Filters gene expression by genes and donors')

    parser.add_argument('--gx',
                         help='input file with gene expression (genes x donors)',
                         type=str,
                         dest='gx_file')

    parser.add_argument('--donors',
                         dest='donor_file',
                         type=str,
                         help='file with donor ids')

    parser.add_argument('--out',
                        dest='out_file',
                        default=None,
                        help='output file for new gene expression')    

    parser.add_argument('--dataset',
                        dest='dataset',
                        help='gtex or cardiogenics')    

    opts = parser.parse_args()
    return opts

def read_samples(donorfile):    
    with open(donorfile, 'r') as samfile:
        sample = 0
        samplenames = list()
        # skip first two lines
        next(samfile)
        next(samfile)
        for line in samfile:
            if re.search('^#', line):
                continue
            sample += 1
            samplenames.append(line.strip().split()[0])
    return samplenames

def filter_donors(df, donors):
    donor_list = df.columns
    common  = [x for x in donors if x in donor_list]
    print("{:d} donors remained from {:d}".format(len(common), len(donor_list)))
    return df[common]

def filter_rows(df, genedict):
    gx_gene_list = df.index
    common  = [genedict[x] for x in gx_gene_list]
    print("{:d} genes remained from {:d}".format(sum(common), len(gx_gene_list)))
    return df[common]

if __name__ == '__main__':
    args = parse_args()

    print("Gene Expr File: {:s}".format(args.gx_file))

    gtfpath = "/cbscratch/franco/datasets/gtex/gencode.v19.annotation.gtf.gz"
    # can't read this with current library
    # gtfpath = "/cbscratch/franco/datasets/gtex/gencode.v28lift37.annotation.gtf.gz"
    print("Reading GENCODE file")

    donors = read_samples(args.donor_file)

    if args.dataset == "gtex":
        gene_info = readgtf.gencode_v12(gtfpath, trim=False)
    if args.dataset == "cardiogenics":
        gene_info = readgtf.gencode_v12(gtfpath, trim=True)

    gene_dict = defaultdict(lambda: False)
    for g in gene_info:
        gene_dict[g.ensembl_id] = True 	

    gx_df = pd.read_table(args.gx_file, sep="\t", header=0, index_col=0)
    new_gx_df = filter_rows(gx_df, gene_dict)
    sorted_gx_df = filter_donors(new_gx_df, donors)
    # outfile = GXFILE+".gencode_filter"
    if args.out_file is None:
        args.out_file = args.gx_file + ".gencode_filtered"
    sorted_gx_df.to_csv(args.out_file, doublequote=False, sep="\t")