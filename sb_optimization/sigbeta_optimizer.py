import numpy as np
import collections
import gzip
import argparse
from scipy import stats
import revreg


SNPINFO_FIELDS = ['chrom', 'varid', 'bp_pos', 'ref_allele', 'alt_allele', 'maf']
class SnpInfo(collections.namedtuple('_SnpInfo', SNPINFO_FIELDS)):
    __slots__ = ()


SEARCHINFO_FIELDS = ['sbeta', 'ngp', 'qvarstd']
class SearchInfo(collections.namedtuple('_SearchInfo', SEARCHINFO_FIELDS)):
    __slots__ = ()


def parse_args():

    parser = argparse.ArgumentParser(description="Given a gene expression file, determine the best sigma beta")

    parser.add_argument('--gx',
                        type=str,
                        dest='gxfile',
                        help='path of input gene expression file')

    parser.add_argument('--out',
                        type=str,
                        dest='outfile',
                        help='path of output file')

    parser.add_argument('--nsnp',
                        type=int,
                        dest='nsnp',
                        default=5000,
                        help='number of SNPs to simulate')
    
    parser.add_argument('--samples',
                        type=str,
                        dest='samplefile',
                        default=None,
                        help='set of samples to use')

    opts = parser.parse_args()
    return opts


def create_snps(nsample, nsnp = 5000, fmin = 0.01, fmax = 0.1):

    f = np.random.uniform(fmin, fmax, nsnp)
    dosage = np.zeros((nsnp, nsample))
    snpinfo = list()
    
    for i in range(nsnp):
        mafratios = np.array([(1 - f[i])**2, 2 * f[i] * (1 - f[i]), f[i]**2])
        accept = False
        while not accept:
            nfreq  = np.random.multinomial(nsample, mafratios, size=1)[0]
            if np.sum(nfreq != 0) >= 2:
                accept = True
        f1 = np.repeat(0, nfreq[0])
        f2 = np.repeat(1, nfreq[1])
        f3 = np.repeat(2, nfreq[2])
        x  = np.concatenate((f1,f2,f3))
        dosage[i, :] = np.random.permutation(x)
        this_snp = SnpInfo(chrom      = 1,
                           bp_pos     = i,
                           varid      = 'rs{:d}'.format(i+1),
                           ref_allele = 'A',
                           alt_allele = 'T',
                           maf        = f[i])
        snpinfo.append(this_snp)

    return dosage, snpinfo


def non_gaussian_parameter(x):
    x2 = np.square(x)
    x4 = np.square(x2)
    data_m4 = np.sum(x4) / x.shape[0]
    stdnorm_m4 = 3
    return (data_m4 / stdnorm_m4) - 1


def sb_optimfunc(sbeta, gt, gx, U, S, Vt):
    pvals, qstat, qmean, qvars, qscale = reverse_regression(gt, gx, sbeta, U, S, Vt)
    kscore = non_gaussian_parameter(qscale)
    vscore = np.std(qvars)
    return kscore, vscore


def reverse_regression(gt_knn, gx_knn, sbeta, U, S, Vt):
    nsnps = gt_knn.shape[0]
    pvals = list()
    qstat = list()
    qmean = list()
    qvars = list()
    qscale = list()

    sigmabeta2 = sbeta * sbeta
    #Yt = gx_knn.T
    #U, S, Vt = np.linalg.svd(Yt, full_matrices=False)
    S2 = np.square(S)

    for i in range(nsnps):
        gt = gt_knn[i, :].copy()
        sigmax2 = np.var(gt)
        S2mod = S2 + sigmax2 / sigmabeta2

        W = np.dot(U, np.dot(np.diag(S2 / S2mod), U.T)) / sigmax2
        Qscore = np.sum(np.square(np.dot(U.T, gt)) * S2 / S2mod) / sigmax2
        pval, muQ, sigmaQ = revreg.pvals_perm(gt.reshape(1, -1), Qscore, W)
        scaledQ = (Qscore - muQ) / sigmaQ

        pvals.append(pval)
        qstat.append(Qscore)
        qmean.append(muQ)
        qvars.append(sigmaQ)
        qscale.append(scaledQ)

    return np.array(pvals), np.array(qstat), np.array(qmean), np.array(qvars), np.array(qscale)


if __name__ == '__main__':

    opts = parse_args()

    # Read the expression file and create N snps for checking deviation from Qscale
    gx_full, gene_list, gx_donors = revreg.read_gtex(opts.gxfile, opts.samplefile)
    gtfull, snp_info = create_snps(gx_full.shape[1], nsnp = opts.nsnp)
    gt_donors = gx_donors
    gx_norm = revreg.normalize_expr(gx_full)
    gx_corr, gt_corr, knn_neighbors = revreg.knn_correction(gx_norm.T, gtfull)
    gx_knn = revreg.normalize_expr(gx_corr.T) #/ np.sqrt(nsample)
    gt_knn = revreg.normalize_and_center_dosage(gt_corr, snp_info)

    Yt = gx_knn.T
    U, S, Vt = np.linalg.svd(Yt, full_matrices=False)

    sbetalist = np.logspace(-3, 0, 100)
    searchpath = list()
    with open(opts.outfile, 'w') as outstream:
        for sbeta in sbetalist:
            kscore, varscore = sb_optimfunc(sbeta, gt_knn, gx_knn, U, S, Vt)
            optpath = SearchInfo(sbeta = sbeta, ngp = kscore, qvarstd = varscore)
            searchpath.append(optpath)
            print(f'{sbeta:g}\t{kscore:g}\t{varscore:g}')
            outstream.write(f'{sbeta}\t{kscore}\t{varscore}\n')
