import time
import sys, os, collections
import numpy as np
import argparse
import mpmath
from operator import attrgetter


def parse_args():

    parser = argparse.ArgumentParser(description='Check DHS enrichment.')

    parser.add_argument('--snpfile',
                        dest='gtfile',
                        metavar='FILE',
                        help='list of significant SNPs, with a space \{:d\} for chromosome number (3 columns: rsid, pos, pval')

    parser.add_argument('--null-snpfile',
                        dest='null_gtfile',
                        metavar='FILE',
                        help='list of significant SNPs from NULL model, with a space \{:d\} for chromosome number (3 columns: rsid, pos, pval')

    parser.add_argument('--dhsfile',
                        dest='dhsfile',
                        metavar='FILE',
                        help='BED annotation file for DHS regions (3 columns: chrm, start, end')

    parser.add_argument('--annotfile',
                        dest='annotfile',
                        metavar='FILE',
                        help='output file for the annotated SNPs in the analysis')

    parser.add_argument('--outfile',
                        dest='outfile',
                        metavar='FILE',
                        help='output file for plots')

    parser.add_argument('--desc',
                        dest='desc',
                        metavar='FILE',
                        type=str,
                        help='description of expression, corrections, etc (eg, tmm_lasso knn10 multi_tissue')

    opts = parser.parse_args()
    return opts

mpmath.mp.dps = 500
def pval(x): return float(mpmath.log10(1 - 0.5 * (1 + mpmath.erf(x/mpmath.sqrt(2)))))

SNPRES_FIELDS = ['rsid', 'chrom', 'pos', 'logp', 'inDHS']
class SNPRes(collections.namedtuple('_SNPRes', SNPRES_FIELDS)):
    __slots__ = ()
        
def tejaas(filepath, chrom, snpannot):
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
            logp  = np.log10(p) if p != 0 else pval( (q - mu) / sigma)
            inDHS = snpannot.get(rsid, False)
            res.append(SNPRes(rsid=rsid, chrom=chrom, pos=pos, logp=-logp, inDHS=inDHS))
    return res

def eqtlgen(filepath, chrom, snpannot):
    res = list()
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr   = line.strip().split("\t")
            rsid  = arr[0]
            pos   = int(arr[1])
            p     = float(arr[2])
            logp  = np.log10(p) if p != 0 else -99
            inDHS = snpannot.get(rsid, False)
            res.append(SNPRes(rsid=rsid, chrom=chrom, pos=pos, logp=-logp, inDHS=inDHS))
    return res

def annotate_snps(dhsfile, gtfile, annot_outfile):
    current_chrm = None
    dhs = open(dhsfile)
    line = dhs.readline()
    f = None
    with open(annot_outfile, 'w') as outf:
        outf.write("chr\trsid\tpos\n")
        while line:
            arr = line.rstrip().split()
            if arr[0][3:] == "X":
                break
            chrm = int(arr[0][3:])
            start = int(arr[1])
            end = int(arr[2])
            if chrm != current_chrm:
                if f:
                    f.close()
                # close previous GT file and open new one
                if current_chrm is not None:
                    tend = time.time()
                    print("CHR{:d} took {:g}s".format(current_chrm, tend-tstart))

                current_chrm = chrm
                if not os.path.exists(gtfile.format(current_chrm)):
                    print("ERR: {:s} does not exist".format(gtfile.format(current_chrm)))
                    continue
                print("Processing CHRM", current_chrm, end=" ")
                f = open(gtfile.format(current_chrm), 'r')
                next(f)
                tstart = time.time()
                gtline = f.readline()
            if not gtline:
                line = dhs.readline()
            while gtline:
                gtarr = gtline.split()
                rsid = gtarr[0]
                pos = int(gtarr[1])
                # print(rsid, pos, start, end)
                if pos < start:
                    gtline = f.readline()
                    continue # go to next snp
                elif pos > end:
                    line = dhs.readline()
                    break # go to next DHS line
                else:
                    # print("-->",chrm, rsid, pos, start, end)
                    outf.write("{:d}\t{:s}\t{:d}\n".format(chrm, rsid, pos))
                    gtline = f.readline()
                    continue
        if f:
            f.close()
            tend = time.time()
            print("CHR{:d} took {:g}s".format(current_chrm, tend-tstart))
        print("Done DHS file")

def read_annot(infile):
    annot_dict = dict()
    with open(infile) as fin:
        _ = fin.readline()
        for line in fin:
            arr = line.rstrip().split()
            annot_dict[arr[1]] = True
    return annot_dict

if __name__ == '__main__':
    
    opts = parse_args()

    gtfile = opts.gtfile
    null_gtfile = opts.null_gtfile
    dhsfile = opts.dhsfile
    annot_outfile = opts.annotfile
    plot_outfile = opts.outfile
    description = opts.desc

    eqtlgen_gtfile = "/cbscratch/franco/datasets/EQTLgen/trans-eQTLs_CHR{:d}"
    eqtlgen_annot  = "/cbscratch/franco/datasets/EQTLgen/SNPs_annots.txt"

    print(gtfile)
    annotate_snps(dhsfile, gtfile, annot_outfile)
    annotate_snps(dhsfile, null_gtfile, annot_outfile+".null")
    annotate_snps(dhsfile, eqtlgen_gtfile, eqtlgen_annot)

    # datadir = "/cbscratch/franco/trans-eqtl/dev-pipeline/gtex_v6_dynamic_knn30/"
    chrms = np.arange(1,23)
    snp_annot = read_annot(annot_outfile)
    snp_annot_null = read_annot(annot_outfile+".null")
    snp_annot_eqtlgen = read_annot(eqtlgen_annot)

    snp_res = list()
    snp_res_null = list()
    snp_res_eqtlgen = list()
    for chrom in chrms:
        print(chrom, end=" ")
        snp_res += tejaas(gtfile.format(chrom), chrom, snp_annot)
        snp_res_null += tejaas(null_gtfile.format(chrom), chrom, snp_annot_null)
        snp_res_eqtlgen += eqtlgen(eqtlgen_gtfile.format(chrom), chrom, snp_annot_eqtlgen)
    print("Done loading SNP results")

    MAX_TOP = 2000
    TOP_W = 100
    snp_list = snp_res
    top_list = sorted(snp_list, key=attrgetter('logp'), reverse=True)
    top_list_eqtlgen = sorted(snp_res_eqtlgen, key=attrgetter('logp'), reverse=True)

    x = np.array([])
    y = np.array([])
    wx = np.array([])
    wy = np.array([])
    ygen = np.array([])
    wygen = np.array([])

    true_eqtls = 0
    rand_eqtls = 0
    eqtlgen_true = 0

    prev = 0
    for TOPN in range(TOP_W, MAX_TOP, TOP_W):
        
        window_true_eqtls = np.sum([x.inDHS for x in top_list[prev:TOPN]])
        true_eqtls += window_true_eqtls

        window_eqtlgen_true = np.sum([x.inDHS for x in top_list_eqtlgen[prev:TOPN]])
        eqtlgen_true += window_eqtlgen_true

        rand_distrib = []
        snp_list_rand = snp_res_null
        top_list_rand = sorted(snp_list_rand, key=attrgetter('logp'), reverse=True)

        window_rand_eqtls = np.sum([x.inDHS for x in top_list_rand[prev:TOPN]])
        rand_eqtls += window_rand_eqtls

        print(true_eqtls, rand_eqtls, true_eqtls/rand_eqtls)
        x = np.append(x, TOPN)
        y = np.append(y, true_eqtls/rand_eqtls)
        ygen = np.append(ygen, eqtlgen_true/rand_eqtls)

        wx = np.append(wx, TOPN)
        wy = np.append(wy, window_true_eqtls/window_rand_eqtls)
        wygen = np.append(wygen, window_eqtlgen_true/window_rand_eqtls)
        prev += TOP_W

    import matplotlib 
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(16,8))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.plot(x, y, label="RR")
    ax1.plot(x, ygen, label="EQTLgen")
    ax1.set_title("Cumulative enrichment")
    ax1.axhline(y=1, color='red')

    ax2.bar(wx, wy, 100, alpha=0.3, label="RR")
    ax2.bar(wx, wygen, 100, alpha=0.3, label="EQTLgen")
    ax2.axhline(y=1, color='red')
    ax2.set_title("Windowed enrichment")

    ax1.set_xlabel("First N SNPS")
    ax2.set_xlabel("First N SNPS")
    ax1.set_ylabel("Enrichment")
    ax2.set_ylabel("Enrichment")

    fig.suptitle(description)
    ax1.legend()
    ax2.legend()

    # fig.tight_layout()
    # fig.subplots_adjust(top=0.88)
    plt.savefig(plot_outfile, bbox_inches='tight')
    print("Plot saved to", plot_outfile)
    plt.show()