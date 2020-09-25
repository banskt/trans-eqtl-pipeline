############################################
# Extracts trans-eQTLs snp genotypes
# and saves it as a dosage file
############################################

import os, sys
sys.path.append("../../")
import collections
import re
import numpy as np
import json
from utils import utils
from utils.readvcf_snp import ReadVCF
import gzip

dosage_outfile = "/cbscratch/franco/datasets/gtex_v8/genotypes/vcfs_SHAPEIT2/non_eQTLs_march2020.dosages.gz"

# gtall_file = "/cbscratch/franco/datasets/gtex_v8/genotypes/gtex_v8.sample"
# def read_samples(samplefile):
#     if os.path.exists(samplefile):
#         with open(samplefile, 'r') as samfile:
#             sample = 0
#             samplenames = list()
#             next(samfile)
#             next(samfile)
#             for line in samfile:
#                 if re.search('^#', line):
#                     continue
#                 sample += 1
#                 samplenames.append(line.strip().split()[0])
#         nsample = sample
#         samplenames = samplenames
#         return samplenames, nsample
    
# allsamples, nall = read_samples(gtall_file)

tissue_file = "../../plots/tissue_table.txt"
json_file   = "../../gtex_v8_metadata.json"
tshorts, tfulls, tstrings = utils.read_tissues_str(tissue_file)
with open(json_file) as instream:
    gtex_meta = json.load(instream)
tissue_colors = dict()
tissue_names = dict()
tissue_nsamples = dict()

for tshort, tfull, tstring in zip(tshorts, tfulls, tstrings):
    if tshort in tshorts:
        tissue_names[tshort] = tstring
        tissue_colors[tshort] = "#" + gtex_meta[tfull]["colorHex"]
        tissue_nsamples[tshort] = gtex_meta[tfull]["rnaSeqSampleCount"]
        

SNPRES_FIELDS = ['rsid']
class SNPRes(collections.namedtuple('_SNPRes', SNPRES_FIELDS)):
    __slots__ = ()

def tejaas(filepath):
    res = list()
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr   = line.strip().split("\t")
            rsid  = arr[0]
            res.append(SNPRes(rsid=rsid))
    return res    
    

basepath = "/usr/users/fsimone/trans-eqtl-pipeline/analysis/jupyter/stability_selection/non_eqtls"
trans_dict = dict()
for tissue in tshorts:
    tejaas_file = os.path.join(basepath, f"{tissue}_non_eqtls.txt")
    
    if not os.path.exists(tejaas_file):
        print("{:s} has no null eqtls".format(tissue))
        continue
    print("Loading ", tissue, end="")
    transeqtls = tejaas(tejaas_file)
    if len(transeqtls) > 0:
        trans_dict[tissue] = transeqtls
        print(" has {:d} non-eqtls".format(len(transeqtls)))
    else:
        trans_dict[tissue] = []
        print(" has 0 non-eqtls")

# Sort trans-eqtls by chromosome to search them in the genotype files

teqtl_varids = list()
for tissue in tshorts:
    teqtl_varids += [snp.rsid for snp in trans_dict[tissue]]
teqtl_varids = list(set(teqtl_varids))

chrm_teqtls = dict()
suma = 0
for chrm in range(1,23):
    chrm_teqtls[chrm] = [x for x in teqtl_varids if x.startswith("chr{:d}_".format(chrm))]
    print(f"chr{chrm} has {len(chrm_teqtls[chrm])} trans-eqtls")
    suma += len(chrm_teqtls[chrm])
print(f"Total trans-eqtls:{suma}")

# Load individual dosage for every trans-eqtl, takes time!

SNPGT_FIELDS = ['varid', 'chrom', 'pos', 'ref', 'alt', 'maf', 'dosage']
class SNPGT(collections.namedtuple('_SNPGT', SNPGT_FIELDS)):
    __slots__ = ()

full_teqtls_gt = collections.defaultdict(dict)
first_donors = list()
for chrm in range(1,23):
    print(f"Reading CHR{chrm}")
    f_vcf = "/cbscratch/franco/datasets/gtex_v8/genotypes/vcfs_SHAPEIT2/0.01/GTEX_v8_2020-02-21_WGS_838Indiv_Freeze.SHAPEIT2_phased_NoMissingGT_SNPfilter_MAF0.01_chr{:d}.vcf.gz".format(chrm)
    samplefile = None
    vcf = ReadVCF(f_vcf, snplist=chrm_teqtls[chrm])
    gtfull = vcf.dosage
    gt_donors = vcf.donor_ids
    if chrm == 1:
        first_donors = gt_donors
    else:
        if not first_donors == gt_donors:
            print("donor error!")
            raise
    snpinfos = vcf.snpinfo
    full_teqtls_gt[chrm] = list()
    for i,snp in enumerate(snpinfos):
        full_teqtls_gt[chrm].append(SNPGT(varid=snp.varid, chrom=snp.chrom, pos=snp.bp_pos, ref=snp.ref_allele, alt=snp.alt_allele, maf=snp.maf, dosage=gtfull[i,:]))

with gzip.open(dosage_outfile, 'wb') as outds:
    for chrm in range(1,23):
        for snp in full_teqtls_gt[chrm]:
            # print(snp.chrom, snp.varid, snp.pos, snp.ref, snp.alt, snp.maf, ' '.join([str(x) for x in snp.dosage]))
            outds.write('{:d} {:s} {:d} {:s} {:s} {:g} {:s}\n'.format(snp.chrom, snp.varid, snp.pos, snp.ref, snp.alt, snp.maf, ' '.join([str(x) for x in snp.dosage])).encode())