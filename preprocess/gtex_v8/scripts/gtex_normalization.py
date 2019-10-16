import pandas as pd
import argparse
import gzip, os
from expression_normalization import qn_normalize, tmm_normalize
from expression_normalization import QC_expression
import numpy as np

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
    if "Description" in df.columns:
        df.drop('Description', axis=1, inplace=True)
    df.index.name = 'gene_id'
    df = df[[i for i in df.columns if '-'.join(i.split('-')[:2]) in donor_ids]]
    return df

def centerscale_expr(Y):
    if isinstance(Y, pd.DataFrame):
        Y_cent = (Y.values - np.mean(Y.values, axis = 1).reshape(-1, 1)) / np.std(Y.values, axis = 1).reshape(-1, 1)
        Y_cent = pd.DataFrame(Y_cent, index=Y.index, columns=Y.columns)
        Y_cent.index.name = Y.index.name
    else:
        Y_cent = (Y - np.mean(Y, axis = 1).reshape(-1, 1)) / np.std(Y, axis = 1).reshape(-1, 1)
    return Y_cent

opts = parse_args()
# expression_threshold=0.1    # 'Selects genes with > expression_threshold expression in at least min_samples')
# count_threshold=5,          # 'Selects genes with > count_threshold reads in at least min_samples')
min_samples=10              # 'Minimum number of samples that must satisfy thresholds')


donor_ids = get_donors(opts.donorspath)
expression_df = read_gct(opts.rpkmpath, donor_ids)
counts_df = read_gct(opts.countspath, donor_ids)

if expression_df.shape[1] < min_samples:
    raise ValueError("tissue has less samples than threshold")

expr_ids = list(expression_df.columns)
tissue_counts_df = counts_df.loc[:,expr_ids]

# match sample ids
newcolumns = ["-".join(i.split("-")[:2]) for i in expr_ids]
expression_df.columns = newcolumns
tissue_counts_df.columns = newcolumns


# QC filtering
# TPM filtering
print('  * QC filtering')
qc_tpm_expr, qc_counts = QC_expression(tissue_counts_df, expression_df)
if not os.path.exists(os.path.join(opts.outdir,"tpms")):
    os.makedirs(os.path.join(opts.outdir,"tpms"))
qc_tpm_expr.to_csv(os.path.join(opts.outdir,"tpms","{:s}_tpms_qcfilter.txt".format(opts.tissue)), sep="\t")
qc_counts.to_csv(os.path.join(opts.outdir,"tpms","{:s}_counts_qcfilter.txt".format(opts.tissue)), sep="\t")

# Apply TMM or QN normalization
print('  * Applying QN')
qn_expr  = centerscale_expr(qn_normalize(qc_tpm_expr))

print('  * Applying TMM')
tmm_expr = centerscale_expr(tmm_normalize(qc_counts))

if not os.path.exists(os.path.join(opts.outdir,"qn")):
    os.makedirs(os.path.join(opts.outdir,"qn"))
if not os.path.exists(os.path.join(opts.outdir,"tmm")):
    os.makedirs(os.path.join(opts.outdir,"tmm"))
qn_expr.to_csv(os.path.join(opts.outdir,"qn","{:s}_qn.txt".format(opts.tissue)), sep="\t")
tmm_expr.to_csv(os.path.join(opts.outdir,"tmm","{:s}_tmm.txt".format(opts.tissue)), sep="\t")
