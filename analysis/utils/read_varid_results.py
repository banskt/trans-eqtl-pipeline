import numpy as np
import collections

EQTL_FIELDS = ['chrom', 'varid', 'bp_pos', 'pval', 'log10pval']
class EqtlInfo(collections.namedtuple('_EQTL_FIELDS', EQTL_FIELDS)):
    __slots__ = ()

def transeqtls(filename, debug = False):
    res = list()
    
    idlist = list()
    chrmlist = list()
    poslist = list()
    pvallist = list()
    
    with open(filename, 'r') as instream:
        for line in instream:
            linesplit = line.strip().split()
            varid = linesplit[0].strip()
            fields = varid.split("_")
            chrm = int(fields[0][3:])
            bppos = int(fields[1])
            pval = 5e-8

            idlist.append(varid)
            chrmlist.append(chrm)
            poslist.append(bppos)
            pvallist.append(pval)

    if len(pvallist) > 0:
        pvalarr = np.array(pvallist)

        ## Sanity check of p-value
        # Are there any nan p-values?
        nan_mask = np.isnan(pvalarr)
        if np.any(nan_mask):
            if debug: print(f'SNPs with nan p-value: {np.sum(nan_mask)}')
            pvalarr[np.where(nan_mask)[0]] = 1.0 ## just ignore this SNP

        # Are there any zero p-values
        zero_mask = pvalarr == 0
        if np.any(zero_mask):
            if debug: print(f'SNPs with zero p-value: {np.sum(zero_mask)}')
            nonzero_pvals = pvalarr[~zero_mask]
            if len(nonzero_pvals) > 0:
                pmin = np.min(nonzero_pvals)
            else:
                pmin = 1e-9
            pvalarr[np.where(zero_mask)[0]] = pmin

        for i, varid in enumerate(idlist):
            this_eqtl = EqtlInfo(chrom = chrmlist[i], varid = varid, bp_pos = poslist[i], 
                                 pval = pvalarr[i], log10pval = -np.log10(pvalarr[i]))
            res.append(this_eqtl)
            
    return res
