#!/usr/bin/env python

import numpy as np
import pandas as pd
import scipy.stats as stats
import warnings
from sklearn import linear_model

def lmcorrect(expression_df, cov_df):   
    donor_ids = expression_df.columns

    #sort donors
    cov_df = cov_df[expression_df.columns]

    reg = linear_model.LinearRegression()
    reg.fit(cov_df.T, expression_df.T)

    # reg.score(df_cov.T, crop_expression_df.T)
    # print(reg.coef_)
    residuals = expression_df - reg.predict(cov_df.T).T
    return residuals, reg.coef_


def normalize_quantiles(df):
    """
    Note: replicates behavior of R function normalize.quantiles from library("preprocessCore")  

    Reference:
     [1] Bolstad et al., Bioinformatics 19(2), pp. 185-193, 2003
    
    Adapted from https://github.com/andrewdyates/quantile_normalize
    """
    
    M = df.values.copy()
    
    Q = M.argsort(axis=0)
    m,n = M.shape

    # compute quantile vector
    quantiles = np.zeros(m)
    for i in range(n):
        quantiles += M[Q[:,i],i]
    quantiles = quantiles / n
    
    for i in range(n):
        # Get equivalence classes; unique values == 0
        dupes = np.zeros(m, dtype=np.int)
        for j in range(m-1):
            if M[Q[j,i],i]==M[Q[j+1,i],i]:
                dupes[j+1] = dupes[j]+1
                
        # Replace column with quantile ranks
        M[Q[:,i],i] = quantiles

        # Average together equivalence classes
        j = m-1
        while j >= 0:
            if dupes[j] == 0:
                j -= 1
            else:
                idxs = Q[j-dupes[j]:j+1,i]
                M[idxs,i] = np.median(M[idxs,i])
                j -= 1 + dupes[j]
        assert j == -1
    
    return M

def edgeR_calcNormFactors(counts_df, ref=None, logratio_trim=0.3, sum_trim=0.05, acutoff=-1e10, verbose=False):
    """
    Calculate TMM (Trimmed Mean of M values) normalization.
    Reproduces edgeR::calcNormFactors.default
    Scaling factors for the library sizes that minimize the log-fold changes between the samples for most genes.
    Effective library size: TMM scaling factor * library size
    References:
     [1] Robinson & Oshlack, 2010
     [2] R functions:
          edgeR::calcNormFactors.default
          edgeR:::.calcFactorWeighted
          edgeR:::.calcFactorQuantile
    """

    # discard genes with all-zero counts
    Y = counts_df.values.copy()
    allzero = np.sum(Y>0,axis=1)==0
    if np.any(allzero):
        Y = Y[~allzero,:]

    # select reference sample
    if ref is None:  # reference sample index
        f75 = np.percentile(Y/np.sum(Y,axis=0), 75, axis=0)
        ref = np.argmin(np.abs(f75-np.mean(f75)))
        if verbose:
            print('Reference sample index: '+str(ref))

    N = np.sum(Y, axis=0)  # total reads in each library

    # with np.errstate(divide='ignore'):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        logR = np.log2((Y/N).T / (Y[:,ref]/N[ref])).T  # log fold change; Mg in [1]
        absE = 0.5*(np.log2(Y/N).T + np.log2(Y[:,ref]/N[ref])).T  # average log relative expression; Ag in [1]
        v = (N-Y)/N/Y
        v = (v.T + v[:,ref]).T  # w in [1]

    ns = Y.shape[1]
    tmm = np.zeros(ns)
    for i in range(ns):
        fin = np.isfinite(logR[:,i]) & np.isfinite(absE[:,i]) & (absE[:,i] > acutoff)
        n = np.sum(fin)

        loL = np.floor(n*logratio_trim)+1
        hiL = n + 1 - loL
        loS = np.floor(n*sum_trim)+1
        hiS = n + 1 - loS
        rankR = stats.rankdata(logR[fin,i])
        rankE = stats.rankdata(absE[fin,i])
        keep = (rankR >= loL) & (rankR <= hiL) & (rankE >= loS) & (rankE <= hiS)
        # in [1], w erroneously defined as 1/v ?
        tmm[i] = 2**(np.nansum(logR[fin,i][keep]/v[fin,i][keep]) / np.nansum(1/v[fin,i][keep]))

    tmm = tmm / np.exp(np.mean(np.log(tmm)))
    return tmm


def edgeR_cpm(counts_df, tmm=None, normalized_lib_sizes=True):
    """
    Return edgeR normalized/rescaled CPM (counts per million)
    Reproduces edgeR::cpm.DGEList
    """
    lib_size = counts_df.sum(axis=0)
    if normalized_lib_sizes:
        if tmm is None:
            tmm = edgeR_calcNormFactors(counts_df)
        lib_size = lib_size * tmm
    return counts_df / lib_size * 1e6
        

# def inverse_quantile_normalization(M):
#     """
#     After quantile normalization of samples, standardize expression of each gene
#     """
#     R = stats.mstats.rankdata(M,axis=1)  # ties are averaged
#     Q = stats.norm.ppf(R/(M.shape[1]+1))
#     return Q

def inverse_normal_transform(M):
    """
    Transform rows to a standard normal distribution
    """
    R = stats.mstats.rankdata(M, axis=1)  # ties are averaged
    if isinstance(M, pd.DataFrame):
        Q = pd.DataFrame(stats.norm.ppf(R/(M.shape[1]+1)), index=M.index, columns=M.columns)
    else:
        Q = stats.norm.ppf(R/(M.shape[1]+1))
    return Q
        
        
def normalize_expression(expression_df, counts_df, expression_threshold=0.1, count_threshold=5, min_samples=10):
    """
    Genes are thresholded based on the following expression rules:
      >=min_samples with >expression_threshold expression values
      >=min_samples with >count_threshold read counts
    """
    # donor_ids = ['-'.join(i.split('-')[:2]) for i in expression_df.columns]
    donor_ids = expression_df.columns
    
    # expression thresholds
    mask = ((np.sum(expression_df>expression_threshold,axis=1)>=min_samples) & (np.sum(counts_df>count_threshold,axis=1)>=min_samples)).values
    
    # apply normalization
    M = normalize_quantiles(expression_df.loc[mask])
    R = inverse_normal_transform(M)

    quant_std_df = pd.DataFrame(data=R, columns=donor_ids, index=expression_df.loc[mask].index)    
    quant_df = pd.DataFrame(data=M, columns=donor_ids, index=expression_df.loc[mask].index)
    return quant_std_df, quant_df

def prepare_expression(rpkm_df, counts_df, expression_threshold = 0.1, count_threshold = 6, min_samples = 30):
    donor_ids = rpkm_df.columns
    mask = ((np.sum(rpkm_df > expression_threshold,axis=1) >= min_samples) & (np.sum(counts_df > count_threshold,axis=1) >= min_samples)).values

    tmm_counts_df = edgeR_cpm(counts_df, normalized_lib_sizes=True)
    norm_df = inverse_normal_transform(tmm_counts_df[mask])

    return norm_df, tmm_counts_df[mask]
