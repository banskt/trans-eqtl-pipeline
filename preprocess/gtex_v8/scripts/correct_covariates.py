import pandas as pd
import argparse
import gzip, os
from expression_normalization import correct_lasso_iterative, lmcorrect
import numpy as np

def parse_args():

    parser = argparse.ArgumentParser(description='Convert rpkm for a given tissue into inverse quantile normalized gene expression')

    parser.add_argument('--input-dir',
                        type=str,
                        dest='indir',
                        metavar='FILE',
                        help='input base directory where all expressions are (parent of qn, tmm, rpkms dirs)')

    parser.add_argument('--tissue',
                        type=str,
                        dest='tissue',
                        metavar='STR',
                        help='exact name of the tissue')

    parser.add_argument('--cov',
                        type=str,
                        dest='covpath',
                        metavar='FILE',
                        help='covariate file')

    parser.add_argument('--outdir',
                        type=str,
                        dest="outdir",
                        metavar='STR',
                        help='output directory')

    opts = parser.parse_args()
    return opts

def centerscale_expr(Y):
    if isinstance(Y, pd.DataFrame):
        Y_cent = (Y.values - np.mean(Y.values, axis = 1).reshape(-1, 1)) / np.std(Y.values, axis = 1).reshape(-1, 1)
        Y_cent = pd.DataFrame(Y_cent, index=Y.index, columns=Y.columns)
        Y_cent.index.name = Y.index.name
    else:
        Y_cent = (Y - np.mean(Y, axis = 1).reshape(-1, 1)) / np.std(Y, axis = 1).reshape(-1, 1)
    return Y_cent

if __name__ == '__main__':
    
    opts = parse_args()

    qn_file = os.path.join(opts.indir,"qn","{:s}_qn.txt".format(opts.tissue))
    tmm_file = os.path.join(opts.indir,"tmm","{:s}_tmm.txt".format(opts.tissue))
    raw_file = os.path.join(opts.indir,"tpms","{:s}_tpms_qcfilter.txt".format(opts.tissue))

    qn_df = pd.read_csv(qn_file, sep="\t", header=0, index_col=0)
    tmm_df = pd.read_csv(tmm_file, sep="\t", header=0, index_col=0)
    raw_df = centerscale_expr(pd.read_csv(raw_file, sep="\t", header=0, index_col=0))

    # load covariates
    covfile = opts.covpath # +"covariates/{:s}_nopeer_covariates_w_age.txt".format(tissue)
    df_cov = pd.read_table(covfile, header=0, index_col=0)

    # scale covariates
    means = np.mean(df_cov.T)
    stds = np.std(df_cov.T)
    diff = (df_cov.T - means) / stds
    scaled_df_cov = diff.T

    qn_lasso, lasso_coefs = correct_lasso_iterative(qn_df, scaled_df_cov)
    tmm_lasso, lasso_coefs = correct_lasso_iterative(tmm_df, scaled_df_cov)
    raw_lasso, lasso_coefs = correct_lasso_iterative(raw_df, scaled_df_cov)

    qn_cclm, coefs = lmcorrect(qn_df, scaled_df_cov)
    tmm_cclm, coefs = lmcorrect(tmm_df, scaled_df_cov)
    raw_cclm, coefs = lmcorrect(raw_df, scaled_df_cov)

    qn_lasso.to_csv(os.path.join(opts.outdir,"qn","{:s}_qn_lasso.txt".format(opts.tissue)), sep="\t")
    tmm_lasso.to_csv(os.path.join(opts.outdir,"tmm","{:s}_tmm_lasso.txt".format(opts.tissue)), sep="\t")
    raw_lasso.to_csv(os.path.join(opts.outdir,"tpms","{:s}_tpms_lasso.txt".format(opts.tissue)), sep="\t")

    qn_cclm.to_csv(os.path.join(opts.outdir,"qn","{:s}_qn_cclm.txt".format(opts.tissue)), sep="\t")
    tmm_cclm.to_csv(os.path.join(opts.outdir,"tmm","{:s}_tmm_cclm.txt".format(opts.tissue)), sep="\t")
    raw_cclm.to_csv(os.path.join(opts.outdir,"tpms","{:s}_tpms_cclm.txt".format(opts.tissue)), sep="\t")