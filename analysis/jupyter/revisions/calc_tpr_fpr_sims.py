import os, sys
import collections
import re
import numpy as np
import json
import pandas as pd
import argparse

# basedir = "/cbscratch/franco/trans-eqtl/simulation/ldpruned2/15064_450_10_800_30_100_100_0.01_0.5_0.0_1.0_0.6_4.0_0.1_20_0.02/"

def parse_args():

    parser = argparse.ArgumentParser(description='Calculate tpr and fpr for simulations')

    parser.add_argument('--inputdir',
                        type=str,
                        dest='inputdir',
                        help='input directory')

    parser.add_argument('--simname',
                        type=str,
                        dest='simname',
                        help='simulation name. eg sim001')

    parser.add_argument('--outdir',
                        type=str,
                        dest='outdir',
                        help='output directory')

    opts = parser.parse_args()
    return opts

def load_transtarget_genes_jpa(simdir, paramsdir):
    target_genes_pvals = collections.defaultdict(lambda: collections.defaultdict(dict))
    target_gene_dir = os.path.join(simdir, paramsdir)
    targfile = os.path.join(target_gene_dir, "all_gene_snp_list.txt")
    if os.path.exists(targfile):
        with open(targfile) as fin:
            next(fin)
            for line in fin:
                arr = line.split("\t")
                variant_id = arr[1]
                ensembl_id = arr[0]
                pval = float(arr[2])
                target_genes_pvals[variant_id][ensembl_id] = pval
    else:
        print("file not found", targfile)
    return target_genes_pvals

def load_transtarget_genes(simdir, paramsdir, mode, cis_dict):
    target_genes_pvals = collections.defaultdict(lambda: collections.defaultdict(dict))
    target_gene_dir = os.path.join(simdir, paramsdir)
    if mode == "knn":
        targfile = os.path.join(target_gene_dir,"gene_snp_list_knn.txt")
    if mode == "cclm":
        targfile = os.path.join(target_gene_dir,"gene_snp_list.txt")
    if os.path.exists(targfile):
        with open(targfile) as fin:
            next(fin)
            for line in fin:
                arr = line.split("\t")
                variant_id = arr[1]
                ensembl_id = arr[0]
                if cis_dict[variant_id] != ensembl_id:
                    pval = float(arr[2])
                    target_genes_pvals[variant_id][ensembl_id] = pval
    return target_genes_pvals

def load_rr_result(simdir, paramsdir):
    trans_eqtls_dict = collections.defaultdict(lambda: collections.defaultdict(dict))
    ff = os.path.join(simdir, paramsdir, "rr.txt")
    if os.path.exists(ff):
        with open(ff) as fin:
            next(fin)
            for line in fin:
                arr = line.split("\t")
                snp = arr[0]
                maf = float(arr[3])
                pval = float(arr[7])
                trans_eqtls_dict[snp]["pval"] = pval
                trans_eqtls_dict[snp]["maf"] = pval
    else:
        print("Error, file does not exist", ff)
        raise
    return trans_eqtls_dict

def load_jpa_result(simdir, paramsdir):
    trans_eqtls_dict = collections.defaultdict(lambda: collections.defaultdict(dict))
    ff = os.path.join(simdir, paramsdir, "all_jpa.txt")
    if os.path.exists(ff):
        with open(ff) as fin:
            next(fin)
            for line in fin:
                arr = line.split("\t")
                snp = arr[0]
                jpa = float(arr[1])
                pval = float(arr[2])
                trans_eqtls_dict[snp]["pval"] = pval
    else:
        print("Error, file does not exist", ff)
        raise
    return trans_eqtls_dict

def load_meqtl_result(simdir, paramsdir):
    trans_eqtls_dict = collections.defaultdict(lambda: collections.defaultdict(dict))
    ff = os.path.join(simdir, paramsdir, "trans_eqtl.txt")
    if os.path.exists(ff):
        with open(ff) as fin:
            next(fin)
            for line in fin:
                arr = line.split("\t")
                snp = arr[0]
                geneid = arr[1]
                pval = float(arr[4])
                fdr  = float(arr[5])
                if snp in trans_eqtls_dict:
                    continue
                trans_eqtls_dict[snp]["pval"] = pval
                trans_eqtls_dict[snp]["fdr"]  = fdr
    else:
        print("Error, file does not exist", ff)
        raise
    return trans_eqtls_dict

def load_meqtl_targetgenes(simdir, paramsdir):
    target_genes_pvals = collections.defaultdict(lambda: collections.defaultdict(dict))
    ff = os.path.join(simdir, paramsdir, "trans_eqtl.txt")
    if os.path.exists(ff):
        with open(ff) as fin:
            next(fin)
            for line in fin:
                arr = line.split("\t")
                snp = arr[0]
                geneid = arr[1]
                pval = float(arr[4])
                fdr  = float(arr[5])
                target_genes_pvals[snp][geneid] = pval
    else:
        print("Error, file does not exist", ff)
        raise
    return target_genes_pvals

# Sorted list of SNPs by pvalue
# sorted([(snp, TranseQTLs_rr["sim001"][snp]["pval"]) for snp in TranseQTLs_rr["sim001"]], key=lambda item: item[1])
def calc_rr_auc(trans_dict_sim, rr_dict):
    tp = 0
    fp = 0
    tpr = list()
    fpr = list()
    sorted_snp_rr = sorted([(snp, rr_dict[snp]["pval"]) for snp in rr_dict], key=lambda item: item[1])
    for snp, pval in sorted_snp_rr:
        if snp in trans_dict_sim:
            tp += 1
        else:
            fp += 1
        tpr.append(tp)
        fpr.append(fp)
    tpr = np.array(tpr)/tp #len(trans_dict_sim.keys())
    fpr = np.array(fpr)/fp #len(rr_dict.keys())
    return tp, fp, tpr, fpr

def calc_meqtl_auc(trans_dict_sim, meqtl_dict):
    tp = 0
    fp = 0
    tpr = list()
    fpr = list()
    sorted_snp = sorted([(snp, meqtl_dict[snp]["pval"]) for snp in meqtl_dict], key=lambda item: item[1])
    for snp, pval in sorted_snp:
        if snp in trans_dict_sim:
            tp += 1
        else:
            fp += 1
        tpr.append(tp)
        fpr.append(fp)
    tpr = np.array(tpr)/len(trans_dict_sim.keys())
    fpr = np.array(fpr)/len(meqtl_dict.keys())
    return tp, fp, tpr, fpr

def read_simtrans(transfile):
    res_dict = collections.defaultdict(dict)
    with open(transfile) as instream:
        for line in instream:
            arr = line.strip().split()
            snpid = arr[1]
            cisgene = arr[3]
            translist = arr[4].split(",")
            transgenes = ["ENSG{:06d}".format(int(x)) for x in translist]
            res_dict[snpid]["targets"] = transgenes
            res_dict[snpid]["cistarget"] = cisgene
    return res_dict

def read_simcis(cisfile):
    res_dict = collections.defaultdict(dict)
    with open(cisfile) as instream:
        for line in instream:
            arr = line.strip().split()
            snpid = arr[1]
            cisgene = arr[3]
            res_dict[snpid] = cisgene
    return res_dict

opts = parse_args()

basedir = opts.inputdir
outdir = opts.outdir
simname = opts.simname
sigmab = 0.2

bTejaas = False
bMatrixeqtl = True
bJPA = True

params_tejaas = f"tejaas/permnull_sb{sigmab}/raw_knn30/peer0"
params_meqtl  = "matrixeqtl/qn_cclm/peer0"
params_jpa    = "jpa/qn_cclm/peer0"

trans_dict = dict()    
simdir = os.path.join(basedir, simname)    
transfile = os.path.join(simdir, "input", "expression.trans")
trans_dict = read_simtrans(transfile)

genepos_file = os.path.join(simdir, "input", "expression.genepos")
genotype_file = os.path.join(simdir, "input", "genotype.dosage.txt")

def map_cis_genes(genepos_file, genotype_file):
    genelist = list()
    snplist  = list()
    cis_dict = dict()
    with open(genepos_file) as fin:
        next(fin)
        for line in fin:
            genelist.append(line.rstrip().split()[0])
    with open(genotype_file) as fin:
        for line in fin:
            snplist.append(line.rstrip().split()[1])
    if len(genelist) == len(snplist):
        cis_dict = dict(zip(snplist, genelist))
    else:
        print("Error, gene and snp numbers are not equal")
        raise
    return cis_dict

cis_dict = map_cis_genes(genepos_file, genotype_file)

TargetedGenes = dict()
TargetedGenes_meqtl = dict()
TargetedGenes_jpa = dict()

TranseQTLs_rr = dict()
TranseQTLs_meqtl = dict()
TranseQTLs_jpa = dict()
modes = ["cclm"] # "knn"


print("Loading trans-eqtl data")
if bTejaas:
    TranseQTLs_rr = load_rr_result(simdir, params_tejaas)
if bMatrixeqtl:
    TranseQTLs_meqtl = load_meqtl_result(simdir, params_meqtl)
if bJPA:
    TranseQTLs_jpa = load_jpa_result(simdir, params_jpa)



print("Loading trans-target data")
if bMatrixeqtl:
    TargetedGenes_meqtl = load_meqtl_targetgenes(simdir, params_meqtl)
if bJPA:
    TargetedGenes_jpa = load_transtarget_genes_jpa(simdir, params_jpa)

for mode in modes:
    if bTejaas:
        TargetedGenes[mode] = load_transtarget_genes(simdir, params_tejaas, mode, cis_dict)
    
print("Running trans-eQTL discovery")

if bTejaas:
    #########################################
    ##### TEJAAS trans-eQTL discovery
    #########################################
    rr_dict = TranseQTLs_rr

    tp, fp, tpr, fpr = calc_rr_auc(trans_dict, rr_dict)
        
    auc_rr = np.trapz(tpr, fpr)
    auc01_rr = np.trapz(tpr[:int(len(tpr)/10)], fpr[:int(len(fpr)/10)])

    with open(os.path.join(outdir, f"rr_tpr_fpr_{simname}_sb{sigmab}.txt"), 'w') as outst:
        outst.write(f"tpr({tp})\tfpr({fp})\n")
        for i in range(len(tpr)):
            outst.write(f"{tpr[i]}\t{fpr[i]}\n")

    with open(os.path.join(outdir, f"rr_auc_{simname}_sb{sigmab}.txt"), 'w') as outst:
        outst.write("auc\tauc01\n")
        outst.write(f"{auc_rr}\t{auc01_rr}\n")

if bJPA:
    #########################################
    ##### TEJAAS trans-eQTL discovery
    #########################################
    jpa_dict = TranseQTLs_jpa

    tp, fp, tpr, fpr = calc_rr_auc(trans_dict, jpa_dict)

    auc_jpa = np.trapz(tpr, fpr)
    auc01_jpa = np.trapz(tpr[:int(len(tpr)/10)], fpr[:int(len(fpr)/10)])  

    with open(os.path.join(outdir, f"jpa_tpr_fpr_{simname}.txt"), 'w') as outst:
        outst.write(f"tpr({tp})\tfpr({fp})\n")
        for i in range(len(tpr)):
            outst.write(f"{tpr[i]}\t{fpr[i]}\n")

    with open(os.path.join(outdir, f"jpa_auc_{simname}.txt"), 'w') as outst:
        outst.write("auc\tauc01\n")
        outst.write(f"{auc_jpa}\t{auc01_jpa}\n")

if bMatrixeqtl:
    #########################################
    ##### MatrixEQTL trans-eQTL discovery
    #########################################
    meqtl_dict = TranseQTLs_meqtl

    tp, fp, tpr, fpr = calc_meqtl_auc(trans_dict, meqtl_dict)
    auc_meqtl = np.trapz(tpr, fpr)
    auc01_meqtl = np.trapz(tpr[:int(len(tpr)/10)], fpr[:int(len(fpr)/10)])

    with open(os.path.join(outdir, f"meqtl_tpr_fpr_{simname}.txt"), 'w') as outst:
        outst.write(f"tpr({tp})\tfpr({fp})\n")
        for i in range(len(tpr)):
            outst.write(f"{tpr[i]}\t{fpr[i]}\n")

    with open(os.path.join(outdir, f"meqtl_auc_{simname}.txt"), 'w') as outst:
        outst.write("auc\tauc01\n")
        outst.write(f"{auc_meqtl}\t{auc01_meqtl}\n")

def bh_procedure_targetgenes(sorted_snp_gene_pval, target_fdr):
    n_tests = len(sorted_snp_gene_pval)
    pass_snps = list()
    for i, snp_pval in enumerate(sorted_snp_gene_pval):
        bh_factor = ((i+1)/n_tests)*target_fdr
        fdr_val = snp_pval[2] / bh_factor
        if snp_pval[2] <= bh_factor:
            pass_snps.append(snp_pval+(fdr_val,))
        else:
            break
    return pass_snps

print("Running trans-target gene discovery")

SIM_FIELDS = ['simname', 'auc', 'auc01', 'auc001', 'cutoff', 'tp', 'tn']
class SimResult(collections.namedtuple('SIMFIELDS', SIM_FIELDS)):
    __slots__ = ()

snp_cutoffs = [0.01, 0.001, 0.0001]


if bTejaas:
    #########################################
    ##### TEJAAS target-gene discovery
    #########################################
    auc_res = list()

    rr_dict = TranseQTLs_rr
    target_dict = TargetedGenes["cclm"]

    ## First, get an ordered list of all snp-gene-pvalues
    snp_gene_pairs = [(snp, geneid, target_dict[snp][geneid]) for snp in target_dict for geneid in target_dict[snp]]
    sorted_pairs = sorted(snp_gene_pairs, key=lambda item: item[2])

    target_fdr = 0.99
    pass_targets = bh_procedure_targetgenes(sorted_pairs, target_fdr)

    pass_targets_checked = list()
    for snp_gene_pval in pass_targets:
        isTeqtl = 0
        isTG = 0
        if snp_gene_pval[0] in trans_dict:
            isTeqtl = 1
            if snp_gene_pval[1] in trans_dict[snp_gene_pval[0]]['targets']:
                isTG = 1
        pass_targets_checked.append(snp_gene_pval + (isTeqtl, isTG))

    with open(os.path.join(outdir, f"rr_target_genes_fdr_bh_{simname}_sb{sigmab}.txt"), 'w') as outst:
        outst.write("snp\tgene\tpval\tfdr\tis_transeqtl\tis_target_gene\n")
        for x in pass_targets_checked:
            outst.write(f"{x[0]}\t{x[1]}\t{x[2]}\t{x[3]}\t{x[4]}\t{x[5]}\n")

    ## then go through that list, but only consider snps where tejaas pval is lower than cutoff
    for snp_cutoff in snp_cutoffs:
        pass_pairs = list()
        for snpid, geneid, pval in sorted_pairs:
            if rr_dict[snpid]['pval'] <= snp_cutoff:
                pass_pairs.append((snpid, geneid, pval))

        print(len(pass_pairs))

        tp = 0
        fp = 0
        tpr = list()
        fpr = list()
        if len(pass_pairs) > 0:
            for snpid, geneid, pval in pass_pairs:
                if snpid in trans_dict:
                    if geneid in trans_dict[snpid]['targets']:
                        tp += 1
                    else:
                        fp += 1
                else:
                    fp += 1
                tpr.append(tp)
                fpr.append(fp)

            print(simname, snp_cutoff, tp, fp)
            tpr_rel = np.array(tpr)/tpr[-1]
            fpr_rel = np.array(fpr)/fpr[-1]
            
            auc = np.trapz(tpr_rel, fpr_rel)
            auc01 = np.trapz(tpr_rel[:int(len(tpr_rel)/10)], fpr_rel[:int(len(fpr_rel)/10)])
            auc001 = np.trapz(tpr_rel[:int(len(tpr_rel)/100)], fpr_rel[:int(len(fpr_rel)/100)])
            res = SimResult(simname=simname, auc=auc, auc01=auc01, auc001=auc001, cutoff=snp_cutoff, tp=tp, tn=fp)

            with open(os.path.join(outdir, f"rr_targetgene_cut{snp_cutoff}_tpr_fpr_{simname}_sb{sigmab}.txt"), 'w') as outst:
                outst.write(f"tpr({tp})\tfpr({fp})\n")
                for i in range(len(tpr_rel)):
                    outst.write(f"{tpr_rel[i]}\t{fpr_rel[i]}\n")
        else:
            res = SimResult(simname=simname, auc=0, auc01=0, auc001=0, cutoff=snp_cutoff, tp=0, tn=0)
        auc_res.append(res)

    with open(os.path.join(outdir, f"rr_targetgene_auc_{simname}_sb{sigmab}.txt"), 'w') as outst:
        outst.write("auc\tauc01\tauc001\tcutoff\tTP\tTN\n")
        for res in auc_res:
            outst.write(f"{res.auc}\t{res.auc01}\t{res.auc001}\t{res.cutoff}\t{res.tp}\t{res.tn}\n")


if bJPA:
    #########################################
    ##### JPA target-gene discovery
    #########################################

    auc_res = list()

    jpa_dict = TranseQTLs_jpa
    target_dict = TargetedGenes_jpa

    ## First, get an ordered list of all snp-gene-pvalues
    snp_gene_pairs = [(snp, geneid, target_dict[snp][geneid]) for snp in target_dict for geneid in target_dict[snp]]
    sorted_pairs = sorted(snp_gene_pairs, key=lambda item: item[2])

    target_fdr = 0.99
    pass_targets = bh_procedure_targetgenes(sorted_pairs, target_fdr)
    pass_targets_checked = list()
    for snp_gene_pval in pass_targets:
        isTeqtl = 0
        isTG = 0
        if snp_gene_pval[0] in trans_dict:
            isTeqtl = 1
            if snp_gene_pval[1] in trans_dict[snp_gene_pval[0]]['targets']:
                isTG = 1
        pass_targets_checked.append(snp_gene_pval + (isTeqtl, isTG))

    with open(os.path.join(outdir, f"jpa_target_genes_fdr_bh_{simname}.txt"), 'w') as outst:
        outst.write("snp\tgene\tpval\tfdr\tis_transeqtl\tis_target_gene\n")
        for x in pass_targets_checked:
            outst.write(f"{x[0]}\t{x[1]}\t{x[2]}\t{x[3]}\t{x[4]}\t{x[5]}\n")

    ## then go through that list, but only consider snps where tejaas pval is lower than cutoff
    for snp_cutoff in snp_cutoffs:
        pass_pairs = list()
        for snpid, geneid, pval in sorted_pairs:
            if jpa_dict[snpid]['pval'] <= snp_cutoff:
                pass_pairs.append((snpid, geneid, pval))

        print(len(pass_pairs))

        tp = 0
        fp = 0
        tpr = list()
        fpr = list()
        if len(pass_pairs) > 0:
            for snpid, geneid, pval in pass_pairs:
                if snpid in trans_dict:
                    if geneid in trans_dict[snpid]['targets']:
                        tp += 1
                    else:
                        fp += 1
                else:
                    fp += 1
                tpr.append(tp)
                fpr.append(fp)

            print(simname, snp_cutoff, tp, fp)
            tpr_rel = np.array(tpr)/tpr[-1]
            fpr_rel = np.array(fpr)/fpr[-1]
            
            auc = np.trapz(tpr_rel, fpr_rel)
            auc01 = np.trapz(tpr_rel[:int(len(tpr_rel)/10)], fpr_rel[:int(len(fpr_rel)/10)])
            auc001 = np.trapz(tpr_rel[:int(len(tpr_rel)/100)], fpr_rel[:int(len(fpr_rel)/100)])
            res = SimResult(simname=simname, auc=auc, auc01=auc01, auc001=auc001, cutoff=snp_cutoff, tp=tp, tn=fp)

            with open(os.path.join(outdir, f"jpa_targetgene_cut{snp_cutoff}_tpr_fpr_{simname}.txt"), 'w') as outst:
                outst.write(f"tpr({tp})\tfpr({fp})\n")
                for i in range(len(tpr_rel)):
                    outst.write(f"{tpr_rel[i]}\t{fpr_rel[i]}\n")
        else:
            res = SimResult(simname=simname, auc=0, auc01=0, auc001=0, cutoff=snp_cutoff, tp=0, tn=0)
        auc_res.append(res)

    with open(os.path.join(outdir, f"jpa_targetgene_auc_{simname}.txt"), 'w') as outst:
        outst.write("auc\tauc01\tauc001\tcutoff\tTP\tTN\n")
        for res in auc_res:
            outst.write(f"{res.auc}\t{res.auc01}\t{res.auc001}\t{res.cutoff}\t{res.tp}\t{res.tn}\n")



if bMatrixeqtl:
    #########################################
    ##### MatrixEQTL target-gene discovery
    #########################################

    meqtl_dict = TranseQTLs_meqtl
    target_dict = TargetedGenes_meqtl

    ## First, get an ordered list of all snp-gene-pvalues
    snp_gene_pairs = [(snp, geneid, target_dict[snp][geneid]) for snp in target_dict for geneid in target_dict[snp]]
    sorted_pairs = sorted(snp_gene_pairs, key=lambda item: item[2])
        
    tp = 0
    fp = 0
    tpr = list()
    fpr = list()
    for snpid, geneid, pval in sorted_pairs:
        if snpid in trans_dict:
            if geneid in trans_dict[snpid]['targets']:
                tp += 1
            else:
                fp += 1
        else:
            fp += 1
        tpr.append(tp)
        fpr.append(fp)

    print(simname, tp, fp)
    tpr_rel = np.array(tpr)/tpr[-1]
    fpr_rel = np.array(fpr)/fpr[-1]

    auc = np.trapz(tpr_rel, fpr_rel)
    auc01 = np.trapz(tpr_rel[:int(len(tpr_rel)/10)], fpr_rel[:int(len(fpr_rel)/10)])
    auc001 = np.trapz(tpr_rel[:int(len(tpr_rel)/100)], fpr_rel[:int(len(fpr_rel)/100)])
    res = SimResult(simname=simname, auc=auc, auc01=auc01, auc001=auc001, cutoff=0, tp=tp, tn=fp)

    with open(os.path.join(outdir, f"meqtl_targetgene_tpr_fpr_{simname}.txt"), 'w') as outst:
        outst.write(f"tpr({tp})\tfpr({fp})\n")
        for i in range(len(tpr_rel)):
            outst.write(f"{tpr_rel[i]}\t{fpr_rel[i]}\n")

    with open(os.path.join(outdir, f"meqtl_targetgene_auc_{simname}.txt"), 'w') as outst:
        outst.write("auc\tauc01\tauc001\tcutoff\tTP\tTN\n")
        outst.write(f"{res.auc}\t{res.auc01}\t{res.auc001}\t{res.cutoff}\t{res.tp}\t{res.tn}\n")