import numpy as np
import collections
import os
import argparse

import sys
sys.path.append('./')
from utils import utils
from utils import read_tejaas_results
from utils import read_eqtlgen_results


def parse_args():

    parser = argparse.ArgumentParser(description='Calculate DHS enrichment for trans-eQTLs')

    parser.add_argument('--dhs',
                        type=str,
                        dest='dhsfile',
                        metavar='FILE',
                        help='full path of the DHS file')

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

    parser.add_argument('--teqtl',
                        type=str,
                        dest='teqtlfile',
                        metavar='FILENAME',
                        help='name of the trans-eQTL file to be analysed')

    parser.add_argument('--eqtlgen',
                        dest='is_eqtlgen',
                        action='store_true',
                        help='Read eQTLGen')

    opts = parser.parse_args()
    return opts


def find_dhs_overlap(res_dict, dhs_file):
    dhs = open(dhs_file)
    line = dhs.readline()
    prev_chrm = 0
    nannot = 0
    while line:
        arr = line.rstrip().split()
        chrm = int(arr[0][3:])
        start = int(arr[1])
        end = int(arr[2])
        if chrm != prev_chrm:
            remaining = res_dict[chrm]
            checked = 0
        if len(remaining) == 0:
            ## No more SNPs in this chromosome, just continue reading the DHS file
            line = dhs.readline()
        else:
            for pos in remaining:
                if pos < start:
                    checked += 1
                    remaining = res_dict[chrm][checked:]
                    continue # go to next SNP
                elif pos > end:
                    line = dhs.readline()
                    break # go to next DHS line, keep checking the remaining results
                else:
                    # this is an annotated SNP
                    checked += 1
                    remaining = res_dict[chrm][checked:]
                    nannot += 1
                    continue # go to next SNP
        prev_chrm = chrm
    dhs.close()
    return nannot


def get_random_dhs_overlap_instance(iterdir, dhs_file):
    res_dict = dict()
    ntot = 0
    for chrm in range(1, 23):
        chrm_bppos_list = list()
        filename = os.path.join(iterdir, f'chr{chrm}.txt')
        with open(filename, 'r') as infile:
            for line in infile:
                line = line.strip()
                spos = int(line.split()[0].strip())
                chrm_bppos_list.append(spos)
                ntot += 1
        chrm_bppos_list.sort()
        res_dict[chrm] = chrm_bppos_list
    n_annotated = find_dhs_overlap(res_dict, dhs_file)
    dhs_frac = float(n_annotated) / float(ntot)
    return dhs_frac


def get_random_dhs_overlap(dirname, dhs_file):
    dhs_frac_rand = list()
    for i in range(10):
        iterdir = os.path.join(dirname, f'random_50000_{i+1 :02d}')
        dhs_frac_rand_instance = get_random_dhs_overlap_instance(iterdir, dhs_file)
        dhs_frac_rand.append(dhs_frac_rand_instance)
        print (f'Iteration {i+1}. Fraction of DHS overlap = {dhs_frac_rand_instance :7.4f}')
    return np.mean(dhs_frac_rand)


def generate_empirical_dhs_enrichment(dirname, dhs_file, dhs_frac_rand):
    enrichment = list()
    for i in range(1000):
        if i % 100 == 0: print (f'Iteration {i}')
        iterdir = os.path.join(dirname, f'random_1000_{i+1 :04d}')
        dhs_frac = get_random_dhs_overlap_instance(iterdir, dhs_file)
        enrichment.append(dhs_frac / dhs_frac_rand)
    return np.array(enrichment)


random_snp_dir = "/usr/users/sbanerj/gtex_v8/genotype/all_samples/random_sampling_SHAPEIT2"
chrmlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

opts = parse_args()
tissuelist = opts.tissuelist
outfile = opts.outfile
dhsfile = opts.dhsfile
resdir = opts.resdir
teqtlfile = opts.teqtlfile
is_eqtlgen = opts.is_eqtlgen

#dhs_frac_rand = get_random_dhs_overlap(random_snp_dir, dhsfile)
#print(f'Fraction of DHS overlap for randomly selected SNPs: {dhs_frac_rand :7.4f}')
#enrichment_rand = generate_empirical_dhs_enrichment(random_snp_dir, dhsfile, dhs_frac_rand)

import pickle
#filehandler = open('dhs_frac_rand_shapeit.pkl', 'wb')
#pickle.dump(dhs_frac_rand, filehandler)
#filehandler = open('enrichment_rand_shapeit.pkl', 'wb')
#pickle.dump(enrichment_rand, filehandler)
dhs_frac_rand = pickle.load( open( "dhs_frac_rand_shapeit.pkl", "rb" ) )
enrichment_rand = pickle.load( open( "enrichment_rand_shapeit.pkl", "rb" ) )
print(f'Fraction of DHS overlap for randomly selected SNPs: {dhs_frac_rand :7.4f}')

fout = open(outfile, 'w')
fout.write(f'TISSUE\tN_TRANSEQTLS\tDHS_FRAC\tENRICHMENT\tP_VALUE\n')

for tissue in tissuelist:

    resfilename = os.path.join(resdir, tissue, teqtlfile)
    if is_eqtlgen:
        transeqtls = read_eqtlgen_results.transeqtls(resfilename)
    else:
        transeqtls = read_tejaas_results.transeqtls(resfilename)
    print(f'{tissue}: {len(transeqtls)} trans-eQTLs')
    nteqtl = len(transeqtls)

    if nteqtl > 0:

        res_dict = dict()
        for chrm in range(1, 23):
            res_dict[chrm] = list()

        for teqtl in transeqtls:
            chrm = teqtl.chrom
            spos = teqtl.bp_pos
            res_dict[chrm].append(spos)

        n_annotated = find_dhs_overlap(res_dict, dhsfile)
        dhs_frac_tissue = n_annotated / nteqtl

        if dhs_frac_tissue > 0:
            dhs_enrichment = dhs_frac_tissue / dhs_frac_rand
            dhs_enrichment_pval = (np.sum(enrichment_rand >= dhs_enrichment) + 1) / (enrichment_rand.shape[0] + 1)
        else:
            dhs_enrichment = 0
            dhs_enrichment_pval = 1.0

        fout.write(f'{tissue}\t{nteqtl}\t{dhs_frac_tissue}\t{dhs_enrichment}\t{dhs_enrichment_pval}\n')

        print (f'Fraction of trans-eQTLs in DHS region: {dhs_frac_tissue}')
        print (f'Enrichment: {dhs_enrichment}')
        print (f'P-Value: {dhs_enrichment_pval}')

    else:
        print("No trans-eQTLs found")

fout.close()
