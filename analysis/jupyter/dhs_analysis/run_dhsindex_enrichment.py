import time
import sys
sys.path.append('../../')
import numpy as np
import collections
import os
import json
import mpmath
from utils import utils
import scipy.stats as ss
import argparse

def parse_args():

    parser = argparse.ArgumentParser(description='Run DHS enrichments')

    # parser.add_argument('--dhs-file',
    #                     type=str,
    #                     dest='infile',
    #                     metavar='FILE',
    #                     help='Input read counts')

    parser.add_argument('--overwrite',
                        dest='overwrite',
                        default=False,
                        action='store_true',
                        help='Overwrite output files')

    opts = parser.parse_args()
    return opts


mpmath.mp.dps = 50
def pvalue(x): return float(mpmath.log10(1 - 0.5 * (1 + mpmath.erf(x/mpmath.sqrt(2)))))

SNPRES_FIELDS = ['rsid', 'chrom', 'pos', 'logp', 'maf']
class SNPRes(collections.namedtuple('_SNPRes', SNPRES_FIELDS)):
    __slots__ = ()
    
type_master_list = ["Cancer / epithelial","Cardiac","Digestive","Lymphoid","Musculoskeletal",\
                    "Myeloid / erythroid","Neural","Organ devel. / renal","Placental / trophoblast",\
                    "Primitive / embryonic","Pulmonary devel.", \
                    "Renal / cancer","Stromal A","Stromal B","Tissue invariant","Vascular / endothelial"]    

def tejaas(filepath):
    res = list()
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr   = line.strip().split("\t")
            rsid  = arr[0]
            chrom = int(arr[1])
            pos   = int(arr[2])
            maf   = float(arr[3])
            q     = float(arr[4])
            mu    = float(arr[5])
            sigma = float(arr[6])
            p     = float(arr[7])
            if sigma == 0:
                continue
            logp  = np.log10(p) if p != 0 else pvalue( (q - mu) / sigma)
            res.append(SNPRes(rsid=rsid, chrom=chrom, pos=pos, logp=-logp, maf=maf))
    return res

def tejaas_old(filepath):
    res = list()
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr   = line.strip().split("\t")
            rsid  = arr[0]
            pos   = int(arr[1])
            p     = float(arr[5])
            chrom = int(arr[6])
            maf   = float(arr[7])
            q     = float(arr[2])
            mu    = float(arr[3])
            sigma = float(arr[4])
            if sigma == 0:
                continue
            logp  = np.log10(p) if p != 0 else pvalue( (q - mu) / sigma)
            res.append(SNPRes(rsid=rsid, chrom=chrom, pos=pos, logp=-logp, maf=maf))
    return res

def read_snplist(filename, mafcutoff=0.01):
    rsidlist = list()
    maflist  = list()
    with open(filename) as instream:
        next(instream)
        for line in instream:
            arr  = line.strip().split("\t")
            rsid = arr[0]
            maf  = float(arr[1])
            if maf >= mafcutoff and maf <= (1 - mafcutoff) :
                rsidlist.append(rsid)
                maflist.append(maf)
    return rsidlist, maflist

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
        arr = line.rstrip().split("\t")
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

def annotated_random_all(gwrsids, dhs_file, isannotated=False):
    res_dict = dict()
    for chrm in range(1, 23):
        res_dict[chrm] = list()
    for var_id in gwrsids:
        info = var_id.split('_')
        chrm = int(info[0][3:])
        bppos = int(info[1])
        res_dict[chrm].append(bppos)
    nannot, nannot_type = find_annotated(res_dict, dhs_file, isannotated)
    return nannot, nannot_type

def annotated_random(gwrsids, nchoose, dhs_file, isannotated=False):
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
    nannot, nannot_type = find_annotated(res_dict, dhs_file, isannotated)
    return nannot, nannot_type

def sample_binomial(n, p, NTIMES):
    array_n = list()
    for i in range(NTIMES):
        n_success = np.random.binomial(n, p)
        array_n.append(n_success)
    return array_n

def enrichment_pval(ntrans, DHS_RANDOM_BG, enrichment):
    randtrans = sample_binomial(ntrans, DHS_RANDOM_BG, 10000000)
    num_null = np.array(randtrans) /  ntrans

    null_enrichments = num_null / DHS_RANDOM_BG
    ecdf = ECDF(null_enrichments)
    pval = 1 - ecdf(enrichment)
    return pval

def sample_rand_bg(snpsfilefmt, dhs_file, nchoose = 20000, niter = 20, isannotated=False):
    if not os.path.exists(snpsfilefmt):
        print("File does not exist", snpsfilefmt)
        raise
    snps_list, maf_list = read_snplist(snpsfilefmt)
    nannot_rand = list()
    nannot_type_list = list()
    nannot_type_array = list()
    print(f'Iteration', end="")
    for k in range(niter):
        nannot_k, nannot_type_k = annotated_random(snps_list, nchoose, dhs_file, isannotated)
        print(f' {k}', end="")
        #print(f'Iteration {k}: {nannot_k}')
        nannot_rand.append(nannot_k)
        nannot_type_array.append([nannot_type_k[k] if k in nannot_type_k else 0 for k in type_master_list ])
    print("")
    frac_rand = np.mean(nannot_rand) / nchoose
    frac_rand_type = np.mean(np.array(nannot_type_array), axis=0) / nchoose
    return frac_rand, len(snps_list), dict(zip(type_master_list, frac_rand_type))

def full_rand_bg(snpsfilefmt, dhs_file, isannotated=False):
    if not os.path.exists(snpsfilefmt):
        print("File does not exist", snpsfilefmt)
        raise
    snps_list, maf_list = read_snplist(snpsfilefmt)
    nannot_rand, nannot_type_array_tmp = annotated_random_all(snps_list, dhs_file, isannotated)
    nannot_type_array = [nannot_type_array_tmp[k] if k in nannot_type_array_tmp else 0 for k in type_master_list]
    frac_rand = nannot_rand / len(snps_list)
    frac_rand_type = np.array(nannot_type_array) / len(snps_list)
    return frac_rand, len(snps_list), dict(zip(type_master_list, frac_rand_type))

def smart_LD_filter(trans_eqtls, trans_eqtls_ld_regions_file):
    pass_snps = list()
    with open(trans_eqtls_ld_regions_file) as instream:
        for line in instream:
            arr = line.strip("\n").split("\t")
            if len(arr[4]) > 0:
                pass_snps.append(arr[0])
    return [x for x in trans_eqtls if x.rsid in pass_snps]


opts = parse_args()
if opts.overwrite:
    print("Overwrite is ON")

tissue_file = "/usr/users/fsimone/trans-eqtl-pipeline/main/tissues.txt"

expressions = ["raw"]
sbtypes = ["sb"]
sbs     = ["0.1"] #["0.005", "0.007"]
Ks       = ["30"]
extras  = [""]
sumdir = "summary_5e-08"
preprocs = list()
use_LD = True
if use_LD:
    print("Use LD results is ON")

basedir  = "/cbscratch/franco/trans-eqtl/dev-pipeline/gtex_v8_lncRNA_freeze"

# title    = "multi_tissue"
# dhs_file = "/cbscratch/franco/datasets/multi-tissue.master.ntypes.simple.hg19_hglift_hg38_clean_sorted.bed"
# isannotated = False

title    = "DHSindex"
dhs_file = "/cbscratch/franco/datasets/DHSindex/DHSindex.bed"
isannotated = True

master_snp_data_file = "/cbscratch/franco/datasets/gtex_v8_dhs_{:s}_background.txt".format(title)
min_trans_eqtl = 0

chrmlist = np.arange(1,23)

json_file = "../../gtex_v8_metadata.json"
tshorts, tfulls = utils.read_tissues(tissue_file)
with open(json_file) as instream:
    gtex_meta = json.load(instream)
tissue_colors  = dict()
tissue_names   = dict()
tissue_samples = dict()
for tshort, tfull in zip(tshorts, tfulls):
    tissue_names[tshort] = tfull
    tissue_colors[tshort] = "#" + gtex_meta[tfull.replace(" ", "_")]["colorHex"]
    tissue_samples[tshort] = gtex_meta[tfull.replace(" ", "_")]["rnaSeqAndGenotypeSampleCount"]
#sorted_tissues = [x[0] for x in sorted(tissue_samples.items(), key=itemgetter(1))]



for sbtype in sbtypes:
    if sbtype == "sbDynamic":
        sb_iter = keffs
    if sbtype == "sb":
        sb_iter = sbs

    for K in Ks:
        for param in sb_iter:
            sb_proc = sbtype+param
            for extra in extras:
                preproc = "permnull_{:s}_knn{:s}{:s}".format(sb_proc, K, extra)
                preprocs.append(preproc)

for preproc in preprocs:
    for expr in expressions:
        resdir = basedir+"/{:s}/{:s}".format(expr, sumdir)
        pcutoffs = [5e-08, 1e-08, 1e-10, 1e-15] #, 1e-08, 1e-10, 1e-15, 1e-20, 1e-25]# [1e-06, 1e-07, 5e-08, 1e-08, 1e-10, 1e-15, 1e-20, 1e-25]
        # trans_eqtls_file = "trans_eqtls_{:s}.txt".format(file_pcutoff)
        trans_eqtls_file = "trans_eqtls.txt"
        if use_LD:
            trans_eqtls_ld_regions_file = trans_eqtls_file+".ld_regions"
            trans_eqtls_file = trans_eqtls_file+".ld_prune"
        snpsfile = "snps_list.txt"

        dhs_frac_rand = dict()
        dhs_frac_type_rand = dict()
        outdir = os.path.join(resdir, "dhs_enrichments", preproc)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        if os.path.exists(master_snp_data_file):
            print("snp_data file exists")
            dhs_frac_type_rand = collections.defaultdict(dict)
            with open(master_snp_data_file) as ifile:
                next(ifile)
                for line in ifile:
                    arr = line.strip().split("\t")
                    if arr[1] == "all":
                        dhs_frac_rand[arr[0]] = float(arr[3])
                    else:
                        dhs_frac_type_rand[arr[0]][arr[1]] = float(arr[3])
        else:
            print("Creating snp_data file")
            with open(master_snp_data_file, 'w') as ofile:
                ofile.write("tissue\ttissue_name\tn_snps\tdhs_frac_rand\n")
                for tissue in tshorts:
                    snpsfilefmt = f'{resdir}/{tissue}/tejaas/{preproc}/{snpsfile}'
                    st = time.time()
                    dhs_frac_rand[tissue], Nsnps, dhs_frac_type_rand[tissue] = sample_rand_bg(snpsfilefmt, dhs_file, nchoose = 20000, niter = 20, isannotated = isannotated)
                    endt = time.time()
                    print("Took ", endt-st)
                    print (f'Fraction of annotated SNPs: {tissue} - all - {dhs_frac_rand[tissue]:7.4f}')
                    ofile.write(f"{tissue}\tall\t{Nsnps}\t{dhs_frac_rand[tissue]}\n")
                    for dhs_type in type_master_list:
                        print(f"{tissue}\t{dhs_type}\t{Nsnps * dhs_frac_type_rand[tissue][dhs_type]}\t{dhs_frac_type_rand[tissue][dhs_type]}")
                        ofile.write(f"{tissue}\t{dhs_type}\t-\t{dhs_frac_type_rand[tissue][dhs_type]}\n")

        for tissue in tshorts:
            if tissue not in dhs_frac_rand:
                # do the sampling and add to file
                print(f"Background for {tissue} not found. Calculating..")
                for pp in preprocs:
                    snpsfilefmt = f'{resdir}/{tissue}/tejaas/{preproc}/{snpsfile}'
                    if not os.path.exists(snpsfilefmt):
                        continue
                    else:
                        dhs_frac_rand[tissue], Nsnps, dhs_frac_type_rand[tissue] = sample_rand_bg(snpsfilefmt, dhs_file, nchoose = 20000, niter = 20, isannotated = isannotated)
                        with open(master_snp_data_file, 'a') as ofile:
                            ofile.write(f"{tissue}\tall\t{Nsnps}\t{dhs_frac_rand[tissue]}\n")
                            for dhs_type in type_master_list:
                                print(f"{tissue}\t{dhs_type}\t{Nsnps * dhs_frac_type_rand[tissue][dhs_type]}\t{dhs_frac_type_rand[tissue][dhs_type]}")
                                ofile.write(f"{tissue}\t{dhs_type}\t-\t{dhs_frac_type_rand[tissue][dhs_type]}\n")
                        break

        for pcutoff in pcutoffs:
            logcutoff = -np.log10(pcutoff)
            transeqtls = dict()
            dhs_annotated = dict()
            dhs_annotated_type = dict()
            enrichment = collections.defaultdict(dict)
            pval_binom = collections.defaultdict(dict)
            for tissue in tshorts:  
                filefmt = f'{resdir}/{tissue}/tejaas/{preproc}/{trans_eqtls_file}'
                if not os.path.exists(filefmt):
                    print("File does not exist", filefmt)
                    continue
                trans_eqtls = tejaas(filefmt)
                if use_LD:
                    # delete the lonely signif SNPs
                    filefmt_regions = f'{resdir}/{tissue}/tejaas/{preproc}/{trans_eqtls_ld_regions_file}'
                    trans_eqtls = smart_LD_filter(trans_eqtls, filefmt_regions)
                trans_eqtls_signif = [x for x in trans_eqtls if x.logp >= logcutoff]
                transeqtls[tissue] = trans_eqtls_signif
                # print(f'{tissue} @ {pcutoff}: {len(transeqtls[tissue])} trans-eQTLs')

                dhs_annotated[tissue] = 0
                if len(transeqtls[tissue]) > min_trans_eqtl:
                    res_dict = dict()
                    for chrm in range(1, 23):
                        res_dict[chrm] = list()
                    for x in  transeqtls[tissue]:
                        res_dict[x.chrom].append(x.pos)
                    dhs_annotated[tissue], dhs_annotated_type[tissue] = find_annotated(res_dict, dhs_file, isannotated)
                    pval_binom[tissue]["all"] = ss.binom_test(dhs_annotated[tissue], len(transeqtls[tissue]), dhs_frac_rand[tissue], alternative='greater')
                    enrichment[tissue]["all"] = float(dhs_annotated[tissue]) / len(transeqtls[tissue]) / dhs_frac_rand[tissue]
                    for dhs_type in type_master_list:
                        pval_binom[tissue][dhs_type] = ss.binom_test(dhs_annotated_type[tissue][dhs_type], len(transeqtls[tissue]), dhs_frac_type_rand[tissue][dhs_type], alternative='greater')
                        enrichment[tissue][dhs_type] = float(dhs_annotated_type[tissue][dhs_type]) / len(transeqtls[tissue]) / dhs_frac_type_rand[tissue][dhs_type]
                    # e_pval[tissue] = enrichment_pval(len(transeqtls[tissue]), dhs_frac_rand[tissue], enrichment[tissue])
                    print (f"{tissue}: {dhs_annotated[tissue]} annotated out of {len(transeqtls[tissue])}. Enrichment = {enrichment[tissue]['all']} - {pval_binom[tissue]['all']}")


            outfile_dhs = os.path.join(outdir, "dhs_enrichment_{:s}_{:g}.txt".format(title, pcutoff))
            if use_LD:
                outfile_dhs = outfile_dhs + ".ld_prune"
            fileperm = "a+"
            write_header = True
            if opts.overwrite:
                fileperm = "w"
            if os.path.exists(outfile_dhs):
                write_header = False
            with open(outfile_dhs, fileperm) as outstream:
                if write_header:
                    outstream.write("tissue\tdhs_type\ttranseqtls\tinDHS\tEnrichment\tpval_binom\n")
                for tissue in tshorts:
                    if tissue in transeqtls:
                        nteqtl = len(transeqtls[tissue])
                        if nteqtl > 0:
                            outstream.write(f"{tissue}\tall\t{nteqtl}\t{dhs_annotated[tissue]}\t{enrichment[tissue]['all']}\t{pval_binom[tissue]['all']}\n")
                            for dhs_type in type_master_list:
                                outstream.write(f"{tissue}\t{dhs_type}\t{nteqtl}\t{dhs_annotated_type[tissue][dhs_type]}\t{enrichment[tissue][dhs_type]}\t{pval_binom[tissue][dhs_type]}\n")