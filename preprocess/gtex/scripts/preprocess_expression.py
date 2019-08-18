import os
import argparse
import numpy as np
import pandas as pd
import collections
from sklearn import linear_model

import rnaseqnorm
import readgtf

def parse_args():

    parser = argparse.ArgumentParser(description = 'Preprocess gene expression')

    parser.add_argument('--tpm',
                        type = str,
                        required = True,
                        dest = 'tpmfile',
                        metavar = 'FILE',
                        help = 'input gene expression TPM file')

    parser.add_argument('--counts',
                        type = str,
                        required = True,
                        dest = 'countsfile',
                        metavar = 'FILE',
                        help = 'input gene expression counts file')

    parser.add_argument('--vcf_sample_list',
                        type = str,
                        dest = 'vcf_sample_list',
                        metavar = 'FILE',
                        help='File listing sample IDs with VCF')

    parser.add_argument('--out',
                        type = str,
                        required = True,
                        dest = 'outfile',
                        metavar = 'FILE',
                        help = 'output gene expression file')

    parser.add_argument('--cov',
                        type = str,
                        required = True,
                        dest = 'covfile',
                        metavar = 'FILE',
                        help = 'input covariates file')


    parser.add_argument('--methods',
                        nargs = '*',
                        default = ['raw'],
                        dest = 'methods',
                        help = 'which method to apply: raw, std, qn, tmm, cclm, cclasso or any combination of them')

    parser.add_argument('--tpm_threshold', 
                        type=np.double,
                        default=0.1, 
                        help='Selects genes with > expression_threshold expression in at least sample_frac_threshold')

    parser.add_argument('--count_threshold',
                        type=np.int32,
                        default=6, 
                        help='Selects genes with >= count_threshold reads in at least sample_frac_threshold samples')

    parser.add_argument('--sample_frac_threshold', 
                        type=np.double, 
                        default=0.2, 
                        help='Minimum fraction of samples that must satisfy thresholds')

    parser.add_argument('--gtf',
                         type=str,
                         help='GTF file to use',
                         dest='gtf_file')

    parser.add_argument('--biotype',
                        nargs = '*',
                        dest = 'biotype',
                        type = str,
                        help = 'type of genes to select')

    opts = parser.parse_args()
    return opts


def read_gct(gct_file):
    """
    Load GCT as DataFrame
    Returns G x N dataframe
    Drops the description, converts full sample names to sample IDs
    """
    df = pd.read_csv(gct_file, sep='\t', skiprows=2, index_col=0)
    df.drop('Description', axis=1, inplace=True)
    df.index.name = 'gene_id'
    colnames = ['-'.join(x.split('-')[:2]) for x in df.columns]
    df.columns = colnames
    return df


def read_covariates(cov_file, col_ix):
    """
    Load covariates as DataFrams
    Returns C x N covariates (C = number of covariates, N = number of samples)
    Orders samples according to 'col_ix'
    """
    df = pd.read_csv(cov_file, sep='\t', index_col=0)
    df.index.name = 'covname'
    df = df[col_ix]
    
    ## drop covariate with equal entries (e.g. gender in tissues like ovary, testis, etc.)
    #eq_mask = np.array([np.all(np.isclose(X, X[0])) for X in np.array(df)])
    #if np.any(eq_mask):
    #    print( f'Dropping covariates for equal elements: {", ".join(df.loc[eq_mask].index)}' )
    #    df = df.loc[~eq_mask]
              
    ## identify and drop collinear covariates
    C = df.astype(np.float64).T
    Q, R = np.linalg.qr(C - np.mean(C, axis=0))
    collinear_mask = np.abs(np.diag(R)) < np.finfo(np.float64).eps * C.shape[1]
    if np.any(collinear_mask):
        print( f'Dropping collinear covariates: {", ".join(df.loc[collinear_mask].index)}' )
        df = df.loc[~collinear_mask]
    
    return df


def read_samples(sample_file):
    with open(sample_file, 'r') as infile:
        content = infile.readlines()
    samples = [x.strip() for x in content]
    return samples


def center_and_scale(Y):
    """
    requires  G x N matrix, where G is the number of genes, and Y is the number of samples
    """
    newY = (Y - np.mean(Y, axis = 1).reshape(-1, 1)) / np.std(Y, axis = 1).reshape(-1, 1)
    return newY


def center_and_scale_df(M):
    return pd.DataFrame(center_and_scale(np.array(M)), index=M.index, columns=M.columns)


def covcorrlm(M, W):
    """
    inputs: M is pd.DataFrame, W is np.array
            M of shape G x N (G = number of genes, N = number of samples)
            W of shape C x N (C = number of covariates)
    """
    X = np.array(M)
    Wnorm = center_and_scale(np.array(W))
    linreg = linear_model.LinearRegression()
    linreg.fit(Wnorm.T, X.T)
    Xcorr = X - linreg.predict(Wnorm.T).T
    #return Xcorr
    return pd.DataFrame(Xcorr, index=M.index, columns=M.columns)


def covcorrlasso(M, W, alpha = 0.05):
    """
    inputs: M is pd.DataFrame, W is np.array
            M of shape G x N (G = number of genes, N = number of samples)
            W of shape C x N (C = number of covariates)
    """
    X = np.array(M)
    Wnorm = center_and_scale(np.array(W))
    Xcorr = np.zeros(X.shape)
    for i in range(X.shape[0]):
        lassoreg = linear_model.Lasso(alpha = alpha)
        lassoreg.fit(Wnorm.T, X[i, :])
        Xcorr[i] = X[i, :] - lassoreg.predict(Wnorm.T).T
    #return Xcorr
    return pd.DataFrame(Xcorr, index=M.index, columns=M.columns)


def tmm_normalization(counts_df, mask):
    tmm_counts_df = rnaseqnorm.edgeR_cpm(counts_df, normalized_lib_sizes = True)
    norm_df = rnaseqnorm.inverse_normal_transform(tmm_counts_df[mask])
    return norm_df


def qn_normalization(tpm_df, mask):
    qn_df = rnaseqnorm.normalize_quantiles(tpm_df.loc[mask])
    norm_df = rnaseqnorm.inverse_normal_transform(qn_df)
    return norm_df


def raw_gx(tpm_df, mask):
    R = tpm_df.loc[mask]
    norm_df = pd.DataFrame(np.array(R), index=R.index, columns=R.columns)
    return norm_df


def filter_genes(df, gene_info):
    gene_list = df.index
    keep = [g.ensembl_id for g in gene_info if g.ensembl_id in gene_list]
    return df.loc[keep]


def pp_options(method, gx_df, mask, cov = None):
    # a pythonic dictionary runs through all functions and creates the full dictionary
    # hence if-else is used.
    if method == 'raw':
        print('   > Raw. Only filtering genes')
        res_df = raw_gx(gx_df, mask)
    if method == 'std':
        print('   > Center and scale')
        res_df = center_and_scale_df(gx_df)
    elif method == 'qn':
        print('   > QN')
        res_df = qn_normalization(gx_df, mask)
    elif method == 'tmm':
        print('   > TMM')
        res_df = tmm_normalization(gx_df, mask)
    elif method == 'cclm':
        print('   > Covariate correction with linear model')
        res_df = covcorrlm(gx_df, cov)
    elif method == 'cclasso':
        print('   > Adaptive covariate correction using LASSO')
        res_df = covcorrlasso(gx_df, cov)
    return res_df


if __name__ == '__main__':

    opts = parse_args()

    # Read gene expression
    print ("Reading gene expression TPM file.")
    tpm_df = read_gct(opts.tpmfile)
    gxdonors = list(tpm_df.columns)
    genes = list(tpm_df.index)

    print ("Reading gene expression Counts file.")
    counts_df = read_gct(opts.countsfile)

    """
    Sample lookup. 
    Only include samples which are also in VCF.
    """
    gtdonors = read_samples(opts.vcf_sample_list)
    common_ix = [x for x in gxdonors if x in gtdonors]
    tpm_df = tpm_df[common_ix]
    counts_df = counts_df[common_ix]

    """
    Expression thresholds.
    Genes are thresholded based on the following expression rules:
      TPM >= tpm_threshold in >= sample_frac_threshold * samples
      read counts >= count_threshold in sample_frac_threshold * samples
    """
    ns = tpm_df.shape[1]
    mask = (
        (np.sum(tpm_df >= opts.tpm_threshold, axis = 1) >= opts.sample_frac_threshold * ns) &
        (np.sum(counts_df >= opts.count_threshold, axis = 1) >= opts.sample_frac_threshold * ns)
    ).values

    # Read covariates
    if any(['cc' in x for x in opts.methods]):
        print ("Reading covariates.")
        cov_df = read_covariates(opts.covfile, tpm_df.columns)
    else:
        cov_df = None

    # Perform the QC preprocessing
    gxpp = collections.defaultdict(lambda: None)

    for method in opts.methods:
        print(f'Applying {method}')
        msteps = method.split('_')
        prevstep = 'none'
        curstep = ''

        gxpp['none'] = tpm_df
        print("  Starting with TPM")
        if method.startswith('tmm'):
            print("  No, wait! Starting with counts here.")
            gxpp['none'] = counts_df

        for mstep in msteps:
            curstep = f'{mstep}' if curstep == '' else f'{curstep}_{mstep}'
            print(f'  Step: {curstep}')
            if gxpp[curstep] is None:
                gxpp[curstep] = pp_options(mstep, gxpp[prevstep], mask, cov = cov_df)
            prevstep = curstep

    # Read GTF file for filtering genes
    print (f'Reading GENCODE file to keep {", ".join(opts.biotype)} genes')
    gene_info = readgtf.gencode_v12(opts.gtf_file, trim=False, biotype=opts.biotype)

    # Filter genes and write output
    print ("Writing all outputs.")
    oprefix = os.path.splitext(opts.outfile)[0]
    osuffix = os.path.splitext(opts.outfile)[1]
    biotype_str = '_'.join(opts.biotype)
    for method in opts.methods:
        ofilename = f'{oprefix}_{method}_{biotype_str}{osuffix}'
        print(f' - {method}: {ofilename}')
        QCdf = filter_genes(gxpp[method], gene_info)
        QCdf.to_csv(ofilename, sep='\t')
