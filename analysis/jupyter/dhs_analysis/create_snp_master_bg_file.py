import sys
sys.path.append('../../')
import numpy as np
import collections
import os
from utils import utils
import scipy.stats as ss
import pandas as pd

def find_annotated(res_dict, dhs_file, isannotated=False):
    dhs = open(dhs_file)
    line = dhs.readline()
    prev_chrm = 0
    nannot = 0
    nannot_type = collections.defaultdict(int)

    sorted_res_dict = dict()
    for chrm in range(1,23):
        sorted_res_dict[chrm] = sorted(res_dict[chrm])

    while line:
        arr = line.rstrip().split()
        if arr[0][3:] == "X" or arr[0][3:] == "Y":
            line = dhs.readline()
            continue
        chrm = int(arr[0][3:])
        start = int(arr[1])
        end = int(arr[2])
        if isannotated:
            atype = arr[9]
        if chrm != prev_chrm:
            remaining = sorted_res_dict[chrm]
            checked = 0
        if len(remaining) == 0:
            ## No more SNPs in this chromosome, just continue reading the DHS file
            line = dhs.readline()
        else:
            for pos in remaining:
                if pos < start:
                    checked += 1
                    remaining = sorted_res_dict[chrm][checked:]
                    continue # go to next SNP
                elif pos > end:
                    line = dhs.readline()
                    break # go to next DHS line, keep checking the remaining results
                else:
                    # this is an annotated SNP
                    checked += 1
                    remaining = sorted_res_dict[chrm][checked:]
                    nannot += 1
                    if isannotated:
                        nannot_type[atype] += 1
                    continue # go to next SNP
        prev_chrm = chrm
    dhs.close()
    return nannot, nannot_type

def annotated_random(gwrsids, nchoose, dhs_file):
    chooseidx = np.sort(np.random.choice(len(gwrsids), nchoose, replace = False))
    res_dict = dict()
    for chrm in range(1, 23):
        res_dict[chrm] = list()
    for idx in chooseidx:
        var_id = gwrsids[idx]
        info = var_id.split('_')
        chrm = int(info[0][3:])
        bppos = int(info[1])
        res_dict[chrm].append(bppos)
    nannot, _nannot_type = find_annotated(res_dict, dhs_file)
    return nannot

def sample_rand_bg(snps_list, dhs_file, nchoose = 20000, niter = 20):
    nannot_rand = list()
    print(f'Iteration', end="")
    for k in range(niter):
        nannot_k = annotated_random(snps_list, nchoose, dhs_file)
        print(f' {k}', end="")
        #print(f'Iteration {k}: {nannot_k}')
        nannot_rand.append(nannot_k)
    print("")
    frac_rand = np.mean(nannot_rand) / nchoose
    return frac_rand

title    = "multi_tissue"
dhs_file = "/cbscratch/franco/datasets/multi-tissue.master.ntypes.simple.hg19_hglift_hg38_clean.bed"

basedir  = "/cbscratch/franco/trans-eqtl/dev-pipeline/gtex_v8_SHAPEIT2" #_freeze"
master_snp_data_file = "/cbscratch/franco/datasets/gtex_v8_dhs_{:s}_background_SHAPEIT2.txt".format(title)
# master_snp_data_file = "/cbscratch/franco/datasets/gtex_v8_dhs_{:s}_background_prev.txt".format(title)

all_snp_ids = pd.read_csv(os.path.join("/cbscratch/franco/from_saikat/gtex_v8_202003/all_variants_pvalues_tejaas.txt"), usecols=[0], header=0, sep="\t")
# all_snp_ids = pd.read_csv(os.path.join("/cbscratch/franco/trans-eqtl/dev-pipeline/gtex_v8_lncRNA_freeze/raw/all_variants_pvalues_tejaas.txt"), usecols=[0], header=0, sep="\t")
snps_list = list(all_snp_ids.values.reshape(-1))

if os.path.exists(master_snp_data_file):
    print("snp_data file exists")
else:
    print("Creating snp_data file")
    with open(master_snp_data_file, 'w') as ofile:
        ofile.write("tissue\ttype\tn_snps\tdhs_frac_rand\n")
        dhs_frac_rand = sample_rand_bg(snps_list, dhs_file, nchoose = 20000, niter = 50)
        print (f'Fraction of annotated SNPs: global - {dhs_frac_rand:7.4f}')
        ofile.write(f"global\tall\t{len(snps_list)}\t{dhs_frac_rand}\n")
