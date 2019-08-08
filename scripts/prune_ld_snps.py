import time
import sys, os, collections
import numpy as np
from operator import attrgetter
import argparse

def parse_args():

    parser = argparse.ArgumentParser(description='Filter SNPs in LD.')

    parser.add_argument('--infiles',
                        nargs='+',
                        dest='infiles',
                        metavar='FILE',
                        help='SNPs pval association list files (RR or mEQTL format)')

    parser.add_argument('--chrm',
                        dest='chrm',
                        type=int,
                        help='CHR number to process')

    parser.add_argument('--ldfile',
                        type=str,
                        dest='ldfile',
                        help='LD map file for that chromosome')


    opts = parser.parse_args()
    return opts


def read_ldfile(ldfile):
    ldict = collections.defaultdict(lambda: False)
    with open(ldfile) as instream:
        next(instream)
        for line in instream:
            arr = line.rstrip().split()
            chrm = int(arr[0])
            pos1 = str(arr[1])
            pos2 = str(arr[2])
            n = int(arr[3])
            r2 = float(arr[4])
            if ldict[pos1]:
                ldict[pos1][pos2] = r2
            else:
                ldict[pos1] = collections.defaultdict(lambda: False)
                ldict[pos1][pos2] = r2
                
            if ldict[pos2]:
                ldict[pos2][pos1] = r2
            else:
                ldict[pos2] = collections.defaultdict(lambda: False)
                ldict[pos2][pos1] = r2
    return ldict

import mpmath
mpmath.mp.dps = 500
def pval(x): return float(1 - 0.5 * (1 + mpmath.erf(x/mpmath.sqrt(2))))
   
SNPRES_FIELDS = ['rsid', 'chrom', 'pos', 'q', 'mu', 'sigma', 'p', 'logp']
class SNPRes(collections.namedtuple('_SNPRes', SNPRES_FIELDS)):
    __slots__ = ()
        
def tejaas(filepath, chrom):
    res = list()
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr   = line.strip().split("\t")
            rsid  = arr[0]
            pos   = int(arr[1])
            p     = float(arr[5])
            q     = float(arr[2])
            mu    = float(arr[3])
            sigma = float(arr[4])
            if sigma == 0:
                continue
            p = p if p != 0 else pval( (q - mu) / sigma)
            logp = np.log10(p) # mpmath.log10
            res.append(SNPRes(rsid=rsid, chrom=chrom, pos=pos, q=q, mu=mu, sigma=sigma, p=p, logp=-logp))
    return res

def tejaas_write(snplist, filepath):
    with open(filepath, 'w') as mfile:
        mfile.write("ID\tPos\tQ\tMu\tSigma\tP\n")
        for snp in snplist:
            fmtstring = "{:s}\t{:d}\t{:g}\t{:g}\t{:g}\t{:g}\n"
            mfile.write(fmtstring.format(snp.rsid, snp.pos, snp.q, snp.mu, snp.sigma, snp.p))

def prune_region(region, myldict):
    start = time.time()
    sorted_region = sorted(region, key=attrgetter('logp'), reverse=True)
    rejected = collections.defaultdict(lambda: False)
    accepted = []
    for snp in sorted_region:
        if rejected[str(snp.pos)]:
            continue
        accepted.append(snp)
        if myldict[str(snp.pos)]:
            for r in myldict[str(snp.pos)].keys():
                rejected[r] = True
    took = time.time() - start
    # print("LD prunning took", took)
    return sorted(accepted, key=attrgetter('pos'), reverse=False)

if __name__ == '__main__':
    # tissues = ["gtex-ms"]

    # expressions = ["lasso", "norm"] 
    # randmethods = ["tejaas_rand_"+str(i) for i in range(11,51)]
    # methods = ["tejaas", "tejaas_rand"] #randmethods

    opts = parse_args()

    snp_files = opts.infiles
    print(snp_files)
    # basedir = "/cbscratch/franco/trans-eqtl/dev-pipeline/gtex_v6_norm/"
    # chroms = np.arange(1,23)
    # isworst = True

    start = time.time()
    print("Loading CHR ", opts.chrm, end="")
    myldict = read_ldfile(opts.ldfile)
    took = time.time() - start
    print(" - {:g} seconds".format(took))
        
    # Do the actual pruning on all datasets and stuff
    pruned_snps = list()
    for infile in snp_files:
        print(infile)
        
        # prune snps
        snplist = tejaas(infile, opts.chrm)
        pruned_snps = prune_region(snplist, myldict)
        
         # write pruned snps
        pruned_outfile = infile+".ld_prune"
        tejaas_write(pruned_snps, pruned_outfile)