#!/usr/bin/env python

import numpy as np
import pandas as pd
import scipy.stats as stats

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
        

def inverse_quantile_normalization(M):
    """
    After quantile normalization of samples, standardize expression of each gene
    """
    R = stats.mstats.rankdata(M,axis=1)  # ties are averaged
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
    R = inverse_quantile_normalization(M)

    quant_std_df = pd.DataFrame(data=R, columns=donor_ids, index=expression_df.loc[mask].index)    
    quant_df = pd.DataFrame(data=M, columns=donor_ids, index=expression_df.loc[mask].index)
    return quant_std_df, quant_df
