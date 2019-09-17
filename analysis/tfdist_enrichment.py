import numpy as np
import collections
import os
import math
import argparse

import sys
sys.path.append('./')
from utils import utils
from utils import read_tejaas_results
from utils import mpl_stylesheet


def parse_args():

    parser = argparse.ArgumentParser(description='Calculate tissue-specific transcription factor enrichment')

    parser.add_argument('--tf',
                        type=str,
                        dest='tffile',
                        metavar='FILE',
                        help='list of transcription factors for the given tissue')

    parser.add_argument('--out',
                        type=str,
                        dest='outfile',
                        metavar='FILE',
                        help='name of the output file')

    parser.add_argument('--resdir',
                        type=str,
                        dest='resdir',
                        metavar='DIR',
                        help='name of the result directory')

    parser.add_argument('--tissue',
                        nargs='*',
                        type=str,
                        dest='tissuelist',
                        metavar='TISSUE',
                        help='short name of all tissues to be analysed')

    parser.add_argument('--cis',
                        type=float,
                        dest='cis_window',
                        default=0.1,
                        help='cis window in megabases')

    opts = parser.parse_args()
    return opts


def find_minimum_distance(spos, starts, ends):
    diffs = zip([spos - x for x in starts], [spos - x for x in ends])
    dists = np.array([math.copysign(min(abs(a), abs(b)), a) if a * b > 0 else 0 for a, b in diffs])
    return dists[np.argmin(np.abs(dists))]

GENEINFO_FIELDS = ['name', 'ensembl_id', 'chrom', 'start', 'end']
class GeneInfo(collections.namedtuple('_GeneInfo', GENEINFO_FIELDS)):
    __slots__ = ()

def read_tflist(tffile):
    tflist = list()
    with open(tffile, 'r') as instream:
        next(instream)
        for line in instream:
            linesplit = line.strip().split()
            ensembl = linesplit[0]
            chrom = int(linesplit[1])
            start = int(linesplit[2])
            end = int(linesplit[3])
            name = linesplit[4]
            tflist.append(GeneInfo(name = name, ensembl_id = ensembl, chrom = chrom, start = start, end = end))
    return tflist

def get_random_mindist(chrmlist, outdir, tflist):
    mindist_rand = list()
    for chrm in chrmlist:
        filename = os.path.join(outdir, f'chr{chrm}.txt')
        tfstarts = [x.start for x in tflist if x.chrom == chrm]
        tfends = [x.end for x in tflist if x.chrom == chrm]
        with open(filename, 'r') as infile:
            for line in infile:
                line = line.strip()
                spos = int(line.split()[0].strip())
                xmin = find_minimum_distance(spos, tfstarts, tfends)
                mindist_rand.append(xmin / 1e6)
    mindist_rand = np.array(mindist_rand)
    return mindist_rand

def get_random_cistf(dirname, chrmlist, tflist):
    tf_frac_rand = 0
    for i in range(10):
        iterdir = os.path.join(dirname, f'random_50000_{i+11 :02d}')
        mindist_rand = get_random_mindist(chrmlist, iterdir, tflist)
        this_tf_frac = np.sum(np.abs(mindist_rand) <= cis_window) / mindist_rand.shape[0]
        print(f'Iteration {i}. Fraction of SNPs with cis TFs = {this_tf_frac}')
        tf_frac_rand += this_tf_frac
    tf_frac_rand /= 10
    return tf_frac_rand

def generate_empirical_cistf_dist(dirname, chrmlist, tflist, tf_frac_rand):
    enrichment_rand = list()
    for i in range(1000):
        if i % 100 == 0:
            print(f'Iteration {i}')
        iterdir = os.path.join(dirname, f'random_1000_{i+1 :04d}')
        mindist_rand = get_random_mindist(chrmlist, iterdir, tflist)
        this_tf_frac = np.sum(np.abs(mindist_rand) <= cis_window) / mindist_rand.shape[0]
        this_tf_enrichment = this_tf_frac / tf_frac_rand
        enrichment_rand.append(this_tf_enrichment)
    
    enrichment_rand = np.array(enrichment_rand)
    return enrichment_rand


random_snp_dir = "/usr/users/sbanerj/gtex_v8/genotype/all_samples/random_sampling"
chrmlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

opts = parse_args()

cis_window = opts.cis_window
tissuelist = opts.tissuelist
outfile = opts.outfile

tflist = read_tflist(opts.tffile)
tf_frac_rand = get_random_cistf(random_snp_dir, chrmlist, tflist)
print(f'Fraction of cis TFs for randomly selected SNPs: {tf_frac_rand :7.4f}')
enrichment_rand = generate_empirical_cistf_dist(random_snp_dir, chrmlist, tflist, tf_frac_rand)

fout = open(outfile, 'w')
fout.write(f'TISSUE\tN_TRANSEQTLS\tCISTF_FRAC\tENRICHMENT\tP_VALUE\n')

for tissue in tissuelist:

    resfilename = os.path.join(opts.resdir, tissue, 'trans_eqtls.txt')
    transeqtls = read_tejaas_results.transeqtls(resfilename)
    print(f'{tissue}: {len(transeqtls)} trans-eQTLs')
    nteqtl = len(transeqtls)

    if nteqtl > 0:
        mindist = list()
        tfstarts = dict()
        tfends = dict()
        for chrm in chrmlist:
            tfstarts[chrm] = [x.start for x in tflist if x.chrom == chrm]
            tfends[chrm] = [x.end for x in tflist if x.chrom == chrm]    
        
        for teqtl in transeqtls:
            chrm = teqtl.chrom
            spos = teqtl.bp_pos
            xmin = find_minimum_distance(spos, tfstarts[chrm], tfends[chrm])
            mindist.append(xmin / 1e6)
        
        tf_frac_tissue = len([x for x in mindist if abs(x) <= cis_window]) / nteqtl
    
        if tf_frac_tissue > 0:
            tf_enrichment = tf_frac_tissue / tf_frac_rand
            tf_enrichment_pval = (np.sum(enrichment_rand >= tf_enrichment) + 1) / (enrichment_rand.shape[0] + 1)
        else:
            tf_enrichment = 0
            tf_enrichment_pval = 1.0
    
        fout.write(f'{tissue}\t{nteqtl}\t{tf_frac_tissue}\t{tf_enrichment}\t{tf_enrichment_pval}\n')
    
        print (f'Fraction of cis TF: {tf_frac_tissue}')
        print (f'Enrichment: {tf_enrichment}')
        print (f'P-Value: {tf_enrichment_pval}')
    
    else:
        print("No trans-eQTLs found")

fout.close()
