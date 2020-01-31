import sys
import os
sys.path.append('/usr/users/fsimone/tejaas')
sys.path.append('/usr/users/fsimone/trans-eqtl-pipeline/analysis')
import numpy as np
from utils import readgtf
from utils import utils
import json
import mpmath
import collections
from iotools.readRPKM import ReadRPKM
from utils.readvcf_snp import ReadVCF
from sklearn.decomposition import PCA
from sklearn import linear_model
import argparse


# for t in `grep -v "#" ../../main/tissues.txt | cut -f 2`; do sbatch -p hh --exclusive -N 1 -t 6-00:00:00 -o lasso_${t}.out -e lasso_${t}.err --wrap="$HOME/opt/miniconda/3/envs/env3.6/bin/python RR_lasso.py --tissue ${t} --outdir /cbscratch/franco/trans-eqtl/dev-pipeline/gtex_v8_lncRNA_freeze/lasso_targets"; done

def parse_args():

    parser = argparse.ArgumentParser(description='Do RR Lasso to get set of target genes for GTEx v8.')

    parser.add_argument('--tissue',
                        dest='tissue',
                        type=str,
                        help='Tissue short name')

    parser.add_argument('--outdir',
                        dest='outdir',
                        type=str,
                        help='Output directory')

    parser.add_argument('--gtf',
                        dest='gtf',
                        type=str,
                        help='GENCODE gtf file')

    parser.add_argument('--basedir',
                        dest='basedir',
                        type=str,
                        help='Directory with summary of trans-eqtls')

    parser.add_argument('--vcf',
                        dest='vcf',
                        type=str,
                        help='VCF file')

    parser.add_argument('--fam',
                        dest='fam',
                        type=str,
                        help='fam file')

    parser.add_argument('--gxfile',
                        dest='gxfile',
                        type=str,
                        help='Gene expression matrix file')

    parser.add_argument('--knn',
                        dest='K',
                        type=int,
                        default=30,
                        help='K-neighbours')

    parser.add_argument('--chrm',
                        dest='chrm',
                        type=int,
                        help='CHRM number')


    opts = parser.parse_args()
    return opts


mpmath.mp.dps = 50
def pvalue(x): return float(mpmath.log10(1 - 0.5 * (1 + mpmath.erf(x/mpmath.sqrt(2)))))

SNPRES_FIELDS = ['rsid', 'chrom', 'pos', 'logp', 'target', 'maf']
class SNPRes(collections.namedtuple('_SNPRes', SNPRES_FIELDS)):
    __slots__ = ()

def knn_correction(expr, dosage, K):
    pca = PCA(n_components=min(expr.shape[0], expr.shape[1]))
    pca.fit(expr) # requires N x G
    expr_pca = pca.transform(expr)

    def gene_distance(a, b):
        return np.linalg.norm(a - b)

    nsample = expr.shape[0]
    distance_matrix = np.zeros((nsample, nsample))
    for i in range(nsample):
        for j in range(i+1, nsample):
            dist = gene_distance(expr_pca[i,:], expr_pca[j,:])
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist

    kneighbor = K
    gx_knn = np.zeros_like(expr)
    gt_knn = np.zeros_like(dosage)
    neighbor_list = list()

    for i in range(nsample):
        neighbors = np.argsort(distance_matrix[i, :])[:kneighbor + 1][1:]
        gx_knn[i, :] = expr[i, :] - np.mean(expr[neighbors, :], axis = 0)
        # noisy_neighbors = np.random.choice(neighbors, size = int(2 * kneighbor / 3), replace = False)
        # noisy_neighbors = np.random.choice(neighbors, size = kneighbor, replace = True )
        noisy_neighbors = neighbors
        gt_knn[:, i] = dosage[:, i] - np.mean(dosage[:, noisy_neighbors], axis = 1)
        neighbor_list.append(neighbors)

    return gx_knn, gt_knn

def select_donors(vcf_donors, expr_donors):
    ''' Make sure that donors are in the same order for both expression and genotype
    '''
    common_donors = [x for x in vcf_donors if x in expr_donors]
    vcfmask = np.array([vcf_donors.index(x) for x in common_donors])
    exprmask = np.array([expr_donors.index(x) for x in common_donors])
    return vcfmask, exprmask

def select_genes(info, names):
    ''' Select genes which would be analyzed. 
        Make sure the indices are not mixed up
    '''
    allowed = [x.ensembl_id for x in info]
    common  = [x for x in names if x in allowed]
    genes = [x for x in info if x.ensembl_id in common]
    indices = [names.index(x.ensembl_id) for x in genes]
    return genes, np.array(indices)

def normalize_and_center_dosage(dosage, snpinfo):
    f = [snp.maf for snp in snpinfo]
    f = np.array(f).reshape(-1, 1)
    gtnorm = (dosage - (2 * f)) / np.sqrt(2 * f * (1 - f))
    gtcent = dosage - np.mean(dosage, axis = 1).reshape(-1, 1)
    return gtnorm, gtcent #rr uses gtcent

def tejaas(filepath, mafcutoff=0.01):
    res = list()
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr   = line.strip().split("\t")
            rsid  = arr[0]
            chrom = int(arr[1])
            pos   = int(arr[2])
            maf   = float(arr[3])
            if maf < mafcutoff or maf > (1-mafcutoff):
                continue
            q     = float(arr[4])
            mu    = float(arr[5])
            sigma = float(arr[6])
            p     = float(arr[7])
            if sigma == 0:
                continue
            logp  = np.log10(p) if p != 0 else pvalue( (q - mu) / sigma)
            res.append(SNPRes(rsid=rsid, chrom=chrom, pos=pos, logp=-logp, maf=maf, target=None))
    return res

def tejaas_targets(filepath):
    res = list()
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr   = line.strip().split("\t")
            geneid  = arr[0]
            snpid   = arr[1]
            res.append((snpid, geneid))
    return res

CISMASK_FIELDS = ['rmv_id', 'apply2']
class CisMask(collections.namedtuple('_CisMask', CISMASK_FIELDS)):
    __slots__ = ()

    @property
    def nsnp(self):
        return len(self.apply2)

    def __repr__(self):
        parent_string = super(CisMask, self).__repr__()
        return '{:s}, nsnp = {:d}'.format(parent_string, self.nsnp)

def get_cismasklist(snpinfo, geneinfo, chrom, window=1e6):
    chr_genes_ix = [[] for ichrm in range(22)] 
    chr_genes = [[] for ichrm in range(22)]
    if chrom is not None:
        chr_genes_ix[chrom - 1] = np.array([i for i, g in enumerate(geneinfo) if g.chrom == chrom])
        chr_genes[chrom - 1] = [geneinfo[ix] for ix in chr_genes_ix[chrom - 1]]
    else:
        for ichrm in range(22):
            chr_genes_ix[ichrm] = np.array([i for i, g in enumerate(geneinfo) if g.chrom == ichrm + 1])
            chr_genes[ichrm] = [geneinfo[ix] for ix in chr_genes_ix[ichrm]]
    genemasks = list()
    iprev = 0
    ichrmprev = 0
    for snp in snpinfo:
        pos = snp.bp_pos
        left = pos - window
        right = pos + window
        ichrm = chrom - 1 if chrom is not None else snp.chrom - 1
        iprev_started = False
        if ichrm != ichrmprev:
            iprev = 0
            ichrmprev = ichrm
        thismask = list()
        for i, g in enumerate(chr_genes[ichrm][iprev:]):
            gstart = g.start
            gend = g.end
            if gstart >= left and gstart <= right:
                # thismask.append(iprev + i)
                thismask.append(chr_genes_ix[ichrm][iprev + i])
                if not iprev_started:
                    new_start_iloc = iprev
                    iprev_started = True
            elif gend >= left and gend <= right:
                # thismask.append(iprev + i)
                thismask.append(chr_genes_ix[ichrm][iprev + i])
                if not iprev_started:
                    new_start_iloc = iprev
                    iprev_started = True
            if gstart > right:
                break
        if len(thismask) > 0:
            #genemasks.append(chr_genes_ix[np.array(thismask)])
            #iprev = thismask[0]
            genemasks.append(np.array(thismask))
            iprev = new_start_iloc
        else:
            genemasks.append(np.array([]))
    return genemasks

def compress_cismasklist(genemasks):
    cismasks = list()
    appendmask = False
    endmask = False
    setprev = False
    snplist = list()
    for i, mask in enumerate(genemasks):
        if not setprev:
            prev_mask = mask
            setprev = True
        if np.all(np.array_equal(mask, prev_mask)):
            snplist.append(i)
        else:
            appendmask = True

        if i == len(genemasks) - 1: endmask = True # no more masks to process

        if appendmask:
            thismask = CisMask(rmv_id = prev_mask, apply2 = snplist)
            cismasks.append(thismask)
            snplist = list([i])
            prev_mask = mask
            if not endmask:
                appendmask = False

        if endmask:
            # if not appendmask:
            #     snplist.append(i)
            thismask = CisMask(rmv_id = mask, apply2 = snplist)
            cismasks.append(thismask)

    return cismasks

def reshape_masked_betas(b, mask, ngenes):
    _b = b.reshape(1, ngenes-len(mask.rmv_id))
    paddedBeta = np.zeros( (1, ngenes) )
    inv_ind = np.delete(np.arange(ngenes), mask.rmv_id)
    paddedBeta[:, inv_ind] = _b
    return paddedBeta.reshape(-1)

def write_target_genes(filehandle, rsid, target_genes):
    line = "{:s}\t{:s}\t{:g}\n"
    for b, g in target_genes:
        filehandle.write(line.format(rsid, g.ensembl_id, b))
        
def RR_Lasso(gx, gt, geneinfo, snp, filehandle):
    lm = linear_model.LassoCV(cv=5)
    lm.fit(gx.T, gt)
    reg_score = lm.score(gx.T, gt)
    print("Score: {:g}".format(reg_score))
    betas_reshape = reshape_masked_betas(lm.coef_, mask, EXPR.shape[0])
    target_gene_inds = np.where(np.abs(betas_reshape) > 0)
    target_genes = [(betas_reshape[j], geneinfo[j]) for j in target_gene_inds[0]]
    write_target_genes(filehandle, snp.varid, target_genes)
    return target_genes


if __name__ == '__main__':
        
    opts = parse_args()

    tissue = opts.tissue

    geneinfo_dict = dict()
    gencode_file = opts.gtf
    geneinfo = readgtf.gencode_v12(gencode_file, biotype = ["protein_coding", "lncRNA"])
    for i,g in enumerate(geneinfo):
        geneinfo_dict[g.ensembl_id] = (i,g)

    basedir = opts.basedir
    # basedirLD = "/cbscratch/franco/trans-eqtl/dev-pipeline/gtex_v8_lncRNA_freeze/summary_LD_5e-08"

    vcf_file = opts.vcf #"/cbscratch/franco/datasets/gtex_v8/genotypes/vcfs_0.01/GTEX_v8_2019-07-29_WGS_838Indiv_Freeze_NoMissingGT_SNPfilter_MAF0.01_chr{:d}.vcf.gz"
    fam_file = opts.fam #"/cbscratch/franco/datasets/gtex_v8/genotypes/gtex_v8.sample"
    gx_file  = opts.gxfile #"/cbscratch/franco/trans-eqtl/new_preprocess_aug2019/gtex_v8/expression/tpms/{:s}_tpms_qcfilter.txt.protein_coding_lncRNA_filtered"
    # gx_file_cclm  = "/cbscratch/franco/trans-eqtl/new_preprocess_aug2019/gtex_v8/expression/tmm/{:s}_tmm_cclm.txt.protein_coding_lncRNA_filtered"


    K = opts.K
    transeqtl_file   = os.path.join(basedir, tissue, "tejaas", "trans_eqtls.txt")
    transeqtls = tejaas(transeqtl_file)

    if len(transeqtls) != 0:
        snplist = [x.rsid for x in transeqtls if x.rsid.startswith("chr"+str(opts.chrm)+"_")]
        if len(snplist) > 0:
            rpkm = ReadRPKM(gx_file.format(tissue), "gtex", npca = 0)
            expression = rpkm.expression
            expr_donors = rpkm.donor_ids
            gene_names = rpkm.gene_names

            if not os.path.exists(opts.outdir): os.makedirs(opts.outdir)

            outfilehandle = open(os.path.join(opts.outdir, "chr{:d}_target_genes_{:s}_lasso.txt".format(opts.chrm, tissue)), 'w') 
            # for chrm in range(1, 23):

            vcf = ReadVCF(vcf_file.format(opts.chrm), snplist=snplist, samplefile=fam_file)
            dosage = vcf.dosage
            gt_donor_ids = vcf.donor_ids
            snpinfo = vcf.snpinfo

            vcfmask, exprmask = select_donors(gt_donor_ids, expr_donors)
            genes, indices = select_genes(geneinfo, gene_names) 

            expr      = expression[:, exprmask]
            gt        = dosage[:, vcfmask]

            gx_corr, gt_corr = knn_correction(expr.T, gt, K=K)
            gx_corr_norm = rpkm._normalize_expr(gx_corr.T)
            gt_corr_norm, gt_corr_cent = normalize_and_center_dosage(gt_corr, snpinfo)
            EXPR  = gx_corr_norm

            cismasklist = get_cismasklist(snpinfo, genes, opts.chrm, window=1e6)
            cismaskcomp = compress_cismasklist(cismasklist)

            # rpkm_cclm = ReadRPKM(gx_file_cclm.format(t), "gtex", npca = 0)
            # expression_cclm = rpkm_cclm.expression
            # expr_cclm = expression_cclm[:, exprmask]
            # gt_norm, gt_cent = normalize_and_center_dosage(gt, snpinfo)

            # # KNN on top of CCLM
            # gx_cclm_corr, gt_cclm_corr = knn_correction(expr_cclm.T, gt, K=30)
            # gx_cclm_corr_norm = rpkm._normalize_expr(gx_cclm_corr.T)
            # gt_cclm_corr_norm, gt_cclm_corr_cent = normalize_and_center_dosage(gt_cclm_corr, snpinfo)

            for mask in cismaskcomp:
                usegenes = np.ones(EXPR.shape[0], dtype=bool)
                if mask.rmv_id.shape[0] > 0: usegenes[mask.rmv_id] = False
                gx_masked      = rpkm._normalize_expr(EXPR[usegenes])
                # gx_masked_cclm = EXPR2[usegenes]
                for i in mask.apply2:
                    GT = gt_corr_cent[i,:]
                    targets = RR_Lasso(gx_masked, GT, genes, snpinfo[i], outfilehandle)      
            outfilehandle.close()
        else:
            print("Tissue {:s} has 0 trans-eqtls in CHRM {:d}".format(tissue, opts.chrm))    
    else:
        print("Tissue {:s} has 0 trans-eqtls".format(tissue))