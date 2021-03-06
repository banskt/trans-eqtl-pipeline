import pandas as pd
import argparse
import gzip, os
from expression_normalization import normalize_expression

def parse_args():

    parser = argparse.ArgumentParser(description='Convert rpkm for a given tissue into inverse quantile normalized gene expression')

    parser.add_argument('--rpkm',
                        type=str,
                        dest='rpkmpath',
                        metavar='FILE',
                        help='input RPKM GCT file for a given tissue')

    parser.add_argument('--counts',
                        type=str,
                        dest='countspath',
                        metavar='FILE',
                        help='input counts GCT file for all tissues')

    parser.add_argument('--tissue',
                        type=str,
                        dest='tissue',
                        metavar='STR',
                        help='exact name of the tissue')

    parser.add_argument('--donors',
                        type=str,
                        dest='donorspath',
                        metavar='STR',
                        help='list of genotyped donors')

    parser.add_argument('--outdir',
                        type=str,
                        dest="outdir",
                        metavar='STR',
                        help='output directory')


    opts = parser.parse_args()
    return opts


def get_donors(path):
    donor_ids = list()
    with open(path, 'r') as instream:
        # skip first two lines
        next(instream)
        next(instream)
        for line in instream:
            donor_ids.append(line.strip().split()[0])
    return donor_ids

def read_gct(gct_file, donor_ids):
    """
    Load GCT as DataFrame
    """    
    df = pd.read_csv(gct_file, sep='\t', skiprows=2, index_col=0)
    df.drop('Description', axis=1, inplace=True)
    df.index.name = 'gene_id'
    df = df[[i for i in df.columns if '-'.join(i.split('-')[:2]) in donor_ids]]
    return df

opts = parse_args()
expression_threshold=0.1    # 'Selects genes with > expression_threshold expression in at least min_samples')
count_threshold=5,          # 'Selects genes with > count_threshold reads in at least min_samples')
min_samples=10              # 'Minimum number of samples that must satisfy thresholds')


donor_ids = get_donors(opts.donorspath)
expression_df = read_gct(opts.rpkmpath, donor_ids)
counts_df = read_gct(opts.countspath, donor_ids)

if expression_df.shape[1] < min_samples:
    raise ValueError("tissue has less samples than threshold")

expr_ids = list(expression_df.columns)
tissue_counts_df = counts_df.loc[:,expr_ids]

print('Normalizing using all genes within %i samples ...' % expression_df.shape[1])
quant_std_df, quant_df = normalize_expression(expression_df, tissue_counts_df,
    expression_threshold=expression_threshold, count_threshold=count_threshold, min_samples=min_samples)

newcolumns = ["-".join(i.split("-")[:2]) for i in quant_std_df.columns]
quant_std_df.columns = newcolumns

# write normalized expression file
print("Writing normalized expression for {:s}".format(opts.tissue))
if not os.path.exists(opts.outdir):
    os.makedirs(opts.outdir)
quant_std_df.to_csv(os.path.join(opts.outdir,"{:s}_normalized.txt".format(opts.tissue)), sep='\t')
