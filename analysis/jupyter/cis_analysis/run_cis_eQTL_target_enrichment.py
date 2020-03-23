import numpy as np
import pandas as pd
import os
import sys
import time
import json
sys.path.append('../../')
sys.path.append('/usr/users/fsimone/tejaas')
from utils import readgtf
from utils import utils
import mpmath
import collections
from operator import attrgetter
import gzip
import scipy.stats as ss

mpmath.mp.dps = 50
def pvalue(x): return float(mpmath.log10(1 - 0.5 * (1 + mpmath.erf(x/mpmath.sqrt(2)))))

SNPRES_FIELDS = ['rsid', 'chrom', 'pos', 'logp', 'target', 'maf']
class SNPRes(collections.namedtuple('_SNPRes', SNPRES_FIELDS)):
    __slots__ = ()
    
CT_FIELDS = ['tissue', 'ncis', 'ntrans', 'ncistrans', 'randtrans', 'enrichment', 'pval', 'binom']
class CisTrans(collections.namedtuple('_CisTrans', CT_FIELDS)):
    __slots__ = ()

def load_snp_maf(filepath, tissue):
    snp_maf_dict = collections.defaultdict(lambda:False)
    with open(filepath) as instream:
        next(instream)
        for line in instream:
            arr = line.strip().split("\t")
            snp_maf_dict[arr[0]] = float(arr[1])
    return snp_maf_dict

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

def read_cis(filepath):
    res = list()
    if not os.path.exists(filepath) or os.stat(filepath).st_size == 0:
        print("File empty or does not exist")
        return res
    with gzip.open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr  = line.decode().strip().split("\t")
            rsid = arr[0]
            if rsid.startswith("chrX"):
                continue
            pos = int(rsid.split("_")[1])
            chrom = int(rsid.split("_")[0][3:])
            gene = arr[1].split(":")[-1].split(".")[0]
            maf  = float(arr[5])
            logp = np.log10(float(arr[6]))
            res.append(SNPRes(rsid=rsid, chrom=chrom, pos=pos, logp=-logp, target=gene, maf=maf))
    return res

# for a set of cis and trans-eQTLs, return the cis-trans ids and the cis-eQTLs with its targets
def cross_ref_cis_trans(trans_ids, cis_eqtls):
    cis_ids = list(set([x.rsid for x in cis_eqtls]))
    
    #Intersection between cis-eqtls (MatrixEQTL) and trans-eqtls (TEJAAS)
    cis_trans_eqtls_ids = list(set.intersection(set(trans_ids), set(cis_ids)))
    
    #set up a dict for fast look up later
    cis_trans_dict = dict()
    for x in cis_trans_eqtls_ids:
        cis_trans_dict[x] = True
    
    # List of cis-trans-eqtls with its target gene
    cis_target_eqtls = [x for x in cis_eqtls if cis_trans_dict.get(x.rsid, False)]

    return cis_trans_eqtls_ids, cis_target_eqtls

def crossref_trans_tejaas(transeqtls, cis_eqtls):
    trans_ids = [x.rsid for x in transeqtls]
    a, b = cross_ref_cis_trans(trans_ids, cis_eqtls)
    return a, b

def sample_background_50000_simple(ciseqtls, randompath):
    randtrans = list()
    chroms    = [str(x) for x in np.arange(1,23)]
    for nid in ["{:03d}".format(x) for x in np.arange(1, 11)]:
        Nrand="50000"
        randomfile = randompath+"random_"+Nrand+"_"+nid

        rand_ids = list()
        for chrm in chroms:
            with open(os.path.join(randomfile, "chr{:s}.txt".format(chrm))) as ins:
                rand_ids += [line.rstrip() for line in ins]

        a, b = cross_ref_cis_trans(rand_ids, ciseqtls)
        randtrans.append( len(a) )
    return np.mean(randtrans)

def sample_rand_bg(ciseqtls, snps_list, nchoose = 20000, niter = 20):
    nchoose = 20000 # len(trans_eqtls_signif) 
    niter = 20
    randtrans = list()
    print(f'Iteration', end="")
    for k in range(niter):
        chooseidx = np.sort(np.random.choice(len(snps_list), nchoose, replace = False))
        rand_ids = [snps_list[i] for i in chooseidx]
        # nannot_k = annotated_random(snps_list, nchoose, dhs_file)
        print(f' {k}', end="")
        a, b = cross_ref_cis_trans(rand_ids, ciseqtls)
        randtrans.append( len(a) )
    print("")
    return np.mean(randtrans)

def sample_binomial(n, p, NTIMES):
    array_n = list()
    for i in range(NTIMES):
        n_success = np.random.binomial(n, p)
        array_n.append(n_success)
    return array_n

def smart_LD_filter(trans_eqtls, trans_eqtls_ld_regions_file):
    pass_snps = list()
    with open(trans_eqtls_ld_regions_file) as instream:
        for line in instream:
            arr = line.strip("\n").split("\t")
            if len(arr[4]) > 0:
                pass_snps.append(arr[0])
    return [x for x in trans_eqtls if x.rsid in pass_snps]

json_file = "../../gtex_v8_metadata.json"
with open(json_file) as instream:
    gtex_meta = json.load(instream)
    
tissue_file = "/usr/users/fsimone/trans-eqtl-pipeline/analysis/plots/tissue_table.txt"
tissues, descriptions = utils.read_tissues(tissue_file)
tissue_names   = dict()
tissue_colors  = dict()
tissue_samples = dict()
for tshort, tfull in zip(tissues, descriptions):
    tissue_names[tshort] = tfull
    tissue_colors[tshort] = "#" + gtex_meta[tfull.replace(" ", "_")]["colorHex"]
    tissue_samples[tshort] = gtex_meta[tfull.replace(" ", "_")]["rnaSeqAndGenotypeSampleCount"]

gene_info = readgtf.gencode_v12("/cbscratch/franco/datasets/GENCODE/gencode.v26.annotation.gtf.gz", trim=True)
gene_info_dict = collections.defaultdict(dict)
for gene in gene_info:
    gene_info_dict[gene.chrom][gene.ensembl_id] = gene.typ

# Filter by allowed snps according to MAF
tejaas_expr = "raw"
K = 30 #not used at the moment
# pcutoff = 5e-8
# use_LD    = True
preproc = "permnull_sb0.1_knn30"
MIN_TRANS = 1
MIN_CIS   = 1
gtexportal_dir = "/cbscratch/franco/datasets/gtex_v8/expression/gtex_portal/eQTLs/GTEx_Analysis_v8_eQTL/"

LDs = [True, False]
sb_optims = [True, False]
pcutoffs = ["5e-8", "1e-10"]

for use_LD in LDs:
    for sb_optim in sb_optims:
        for pcutoff in pcutoffs:

            if sb_optim:
                gamma_suffix = "optim_gamma"
            else:
                gamma_suffix = "gamma01"

            if use_LD:
                trans_file = "trans_eqtls_ldpruned.txt"
                out_suffix = f"pc_lncRNA_{gamma_suffix}_knn30_cut{pcutoff}_ldpruned.txt"
                regions_file = "ld_regions.txt"
            else:
                trans_file = "trans_eqtls.txt"
                out_suffix = f"pc_lncRNA_{gamma_suffix}_knn30_cut{pcutoff}.txt"

            # summary_dir = "{:s}/summary_{:g}".format(tejaas_expr,pcutoff)
            # basepath = "/cbscratch/franco/trans-eqtl/dev-pipeline/gtex_v8_lncRNA_freeze/"
            # baseoutdir = os.path.join(basepath, summary_dir, "GTExPortal_eqtl_analysis", preproc)

            basepath = f"/cbscratch/franco/trans-eqtl/protein_coding_lncRNA_{gamma_suffix}_knn30_cut{pcutoff}"
            baseoutdir = os.path.join(basepath, "GTExPortal_eqtl_analysis")
            if not os.path.exists(baseoutdir): os.makedirs(baseoutdir)

            # read all variants ids
            all_snp_ids = pd.read_csv(os.path.join("/cbscratch/franco/from_saikat/gtex_v8_202003/all_variants_pvalues_tejaas.txt"), usecols=[0], header=0, sep="\t")
            snps_list = list(all_snp_ids.values.reshape(-1))

            GENEINFO_FIELDS = ['name', 'ensembl_id', 'chrom', 'start', 'end', 'typ']
            class GeneInfo(collections.namedtuple('_GeneInfo', GENEINFO_FIELDS)):
                __slots__ = ()
                
            def read_TFannot(infile):
                TF_list = list()
                with open(infile) as instream:
                    next(instream)
                    for line in instream:
                        arr = line.rstrip().split()
                        TF_list.append(GeneInfo(ensembl_id=arr[0], chrom=int(arr[1]), start=int(arr[2]), end=int(arr[3]), name=arr[4], typ="TF"))
                return TF_list

            base_dir = "/cbscratch/franco/datasets"

            #### Load all data dictionaries
            TFs_file = "../../external/TF_annotation.txt"
            if not os.path.exists(TFs_file):
                lambert_file = "../../external/TFs_Lambert_2018.csv"
                TFs = pd.read_csv(lambert_file, header = 0)
                TrueTFs = list(TFs[ TFs["isTF?"] == "Yes" ].ID)
                TF_annot = [x for x in gene_info if x.ensembl_id in TrueTFs]
                ## Generate annotation file for future analysis
                with open(TFs_file, 'w') as outstream:
                    outstream.write("ensembl_id\tchrom\tstart\tend\tname\n")
                    for e in TF_annot:
                        outstream.write("\t".join([e.ensembl_id, str(e.chrom), str(e.start), str(e.end), e.name])+"\n")

            TF_annot = read_TFannot(TFs_file)
            TF_dict = collections.defaultdict(dict)
            for g in TF_annot:
                TF_dict[g.chrom][g.ensembl_id] = "TF"

            # Reformat genetype dict, we can add as many gene annotations as we want here
            alltypes_dict = collections.defaultdict(dict)
            genetypes = []
            for chrm in range(1,23):
                gene_info_dict[chrm]
                for k in gene_info_dict[chrm].keys():
                    genetype = gene_info_dict[chrm][k]
                    if genetype not in alltypes_dict:
                        alltypes_dict[genetype] = collections.defaultdict(lambda:False)
                        genetypes.append(genetype)
                    alltypes_dict[genetype][k] = True
                # Add TF dictionary
                for k in TF_dict[chrm].keys():
                    genetype = "TF"
                    if genetype not in alltypes_dict:
                        alltypes_dict[genetype] = collections.defaultdict(lambda:False)
                        genetypes.append(genetype)
                    alltypes_dict[genetype][k] = True
            ##### End data dict load

            ### Calculate Backgrounds
            cis_bg = dict()
            cis_bg_file = os.path.join(gtexportal_dir, "gtex_background_freqs_ciseqtls_SHAPEIT2.txt")
            if not os.path.exists(cis_bg_file):
                for tissue in tissues:
                    print(tissue, end=" ")
                    signif_cisfile = os.path.join(gtexportal_dir, "{:s}.v8.signif_variant_gene_pairs.txt.gz".format(tissue_names[tissue].replace(" ", "_")))
                    if not os.path.exists(signif_cisfile) or os.stat(signif_cisfile).st_size == 0:
                        print("{:s} has no cis-file in GTEx!".format(tissue_names[tissue]))
                        continue
                    ciseqtls = read_cis(signif_cisfile)
                    cis_bg[tissue] = sample_rand_bg(ciseqtls, snps_list, nchoose = 20000, niter = 20)
                with open(cis_bg_file, 'w') as outstream:
                    outstream.write("tissue\trand_count\tfreq\n")
                    for t in cis_bg:
                        outstream.write(f"{t}\t{cis_bg[t]}\t{cis_bg[t]/20000} \n")
            ### End Bg calculation

            res_dict = dict()
            res_target_dict = dict()

            ### Load Backgrounds
            cis_bg_freq = dict()
            with open(cis_bg_file) as instream:
                next(instream)
                for line in instream:
                    arr = line.strip().split("\t")
                    cis_bg_freq[arr[0]] = float(arr[2])

            for tissue in tissues:
                res_dict = dict()
            res_target_dict = dict()
            
            cis_bg_freq = dict()
            with open(cis_bg_file) as instream:
                next(instream)
                for line in instream:
                    arr = line.strip().split("\t")
                    cis_bg_freq[arr[0]] = float(arr[2])

            ################# Calculate Enrichment #################
            for tissue in tissues:   
                # tejaas_file = os.path.join(basepath, summary_dir, tissue, "tejaas", preproc, "trans_eqtls.txt")
                tejaas_file = os.path.join(basepath, tissue, trans_file)
                if not os.path.exists(tejaas_file):
                    print("File does not exist!", tejaas_file)
                    # print("{:s} has no trans-eqtl results".format(tissue))
                    raise
                    #continue
                transeqtls = tejaas(tejaas_file)
                if use_LD:
                    transeqtls = smart_LD_filter(transeqtls, os.path.join(basepath, tissue, regions_file))
                
                if len(transeqtls) < MIN_TRANS:
                    print("{:s} has less than {:d} trans-eqtls".format(tissue, MIN_TRANS))
                    continue
                
                signif_cisfile = os.path.join(gtexportal_dir, "{:s}.v8.signif_variant_gene_pairs.txt.gz".format(tissue_names[tissue].replace(" ", "_")))
                if not os.path.exists(signif_cisfile) or os.stat(signif_cisfile).st_size == 0:
                    print("{:s} has no cis-file in GTEx!".format(tissue_names[tissue]))
                    continue
                ciseqtls = read_cis(signif_cisfile)
                cis_ids = list(set([x.rsid for x in ciseqtls]))
                
                if len(ciseqtls) < MIN_CIS:
                    print("{:s} has less than {:d} cis-eqtls".format(tissue, MIN_CIS))
                    continue
                
                cis_trans_eqtls_ids, cistrans_target_eqtls = crossref_trans_tejaas(transeqtls, ciseqtls)
                
                FRAC_CISTRANS = len(cis_trans_eqtls_ids) / len(transeqtls)
                FRAC_RANDOM_GWCISTRANS = cis_bg_freq[tissue] # randtrans / 50000 
                
                enrichment = FRAC_CISTRANS / FRAC_RANDOM_GWCISTRANS
            
                ncis = len(cis_ids)
                ntrans = len(transeqtls)
                ncistrans = len(cis_trans_eqtls_ids)
                pval_binom = ss.binom_test(ncistrans, ntrans, FRAC_RANDOM_GWCISTRANS, alternative='greater')
                
                res_dict[tissue] = CisTrans(tissue=tissue, ncis=ncis, ntrans=ntrans, 
                                            ncistrans=ncistrans, randtrans=FRAC_RANDOM_GWCISTRANS,
                                            enrichment=enrichment, pval=1, binom=pval_binom)
                
                res_target_dict[tissue] = cistrans_target_eqtls
                
                print(f"########## Tissue: {tissue} - {ntrans} trans-eqtls - {ncistrans} cis-trans-eqtls #########")
                print(f"{tissue:>20}        Enrichment: {enrichment:>g} - binom {pval_binom:>g}")

            ##################### Save tissue enrichments #################

            outcisfilename = os.path.join(baseoutdir,f"ciseqtl_enrichment_"+out_suffix) 
            targets_outfile = os.path.join(baseoutdir,f"ciseqtl_targets_"+out_suffix) 

            # if use_LD:
            #     outcisfilename = os.path.join(baseoutdir,f"ciseqtl_enrichment_pc_lncRNA_{gamma_suffix}_knn30_cut{pcutoff}_ldpruned.txt")
            #     targets_outfile = os.path.join(baseoutdir,f"ciseqtl_targets_pc_lncRNA_{gamma_suffix}_knn30_cut{pcutoff}_ldpruned.txt")
            # else:
            #     outcisfilename = os.path.join(baseoutdir,f"ciseqtl_enrichment_pc_lncRNA_{gamma_suffix}_knn30_cut{pcutoff}.txt")
            #     targets_outfile = os.path.join(baseoutdir,f"ciseqtl_targets_pc_lncRNA_{gamma_suffix}_knn30_cut{pcutoff}.txt")

            with open(outcisfilename, 'w') as outstream:
                for tissue in tissues:
                    if tissue in res_dict:
                        line = f"{tissue}\t{res_dict[tissue].ncis}\t{res_dict[tissue].ntrans}\t{res_dict[tissue].ncistrans}\t{res_dict[tissue].randtrans}\t{res_dict[tissue].enrichment}\t{res_dict[tissue].binom}\n"
                        outstream.write(line)

            with open(targets_outfile, 'w') as outstream:
                for tissue in tissues:
                    if tissue in res_target_dict:
                        for snp in res_target_dict[tissue]:
                            line = f"{tissue}\t{snp.rsid}\t{snp.logp}\t{snp.target}\t{snp.maf}\n"
                            outstream.write(line)
            #################### End save #################

            #################### Read cis-targets and calculate enrichments #################
            res_target_dict = collections.defaultdict(list)
            with open(targets_outfile) as instream:
                for line in instream:
                    arr = line.strip().split("\t")
                    t    = arr[0]
                    rsid = arr[1]
                    logp = -float(arr[2])
                    targ = arr[3]
                    maf  = float(arr[4])
                    chrom = int(rsid.split("_")[0][3:])
                    pos   = int(rsid.split("_")[1])
                    res_target_dict[t].append( SNPRes(rsid=rsid, chrom=chrom, pos=pos, logp=-logp, maf=maf, target=targ) )

            CisTrans_Type_FIELDS = ['tissue', 'genetype', 'ncistrans', 'ntype', 'frac_cis', 'frac_cistrans', 'enrichment', 'binom']
            class CisTrans_Type(collections.namedtuple('_CisTrans_Type', CisTrans_Type_FIELDS)):
                __slots__ = ()

            enrichment_genetypes = list()
            for tissue in tissues:
                signif_cisfile = os.path.join(gtexportal_dir, "{:s}.v8.signif_variant_gene_pairs.txt.gz".format(tissue_names[tissue].replace(" ", "_")))
                if not os.path.exists(signif_cisfile) or os.stat(signif_cisfile).st_size == 0:
                    print("{:s} has no cis-file in GTEx!".format(tissue_names[tissue]))
                    continue
                ciseqtls = read_cis(signif_cisfile)
                
                TOTAL_CIS = len(ciseqtls)
                cis_counts_dict = collections.defaultdict(int)
                for eqtl in ciseqtls:
                    for genetype in genetypes:
                        if alltypes_dict[genetype][eqtl.target]:
                            cis_counts_dict[genetype] += 1
                
                if tissue in res_target_dict:
                    TOTAL_CIS_TRANS = len(res_target_dict[tissue])
                    cistrans_counts_dict = collections.defaultdict(int)
                    for eqtl in res_target_dict[tissue]:
                        for genetype in genetypes:
                            if alltypes_dict[genetype][eqtl.target]:
                                cistrans_counts_dict[genetype] += 1
                            
                    for genetype in genetypes:
                        if TOTAL_CIS > 0 and TOTAL_CIS_TRANS > 0:
                            frac_cis      = cis_counts_dict[genetype] / TOTAL_CIS
                            frac_cistrans = cistrans_counts_dict[genetype] / TOTAL_CIS_TRANS
                            if frac_cis > 0: 
                                if frac_cistrans > 0:
                                    enrichment = frac_cistrans / frac_cis
                                    pval_binom = ss.binom_test(cistrans_counts_dict[genetype], TOTAL_CIS_TRANS, frac_cis, alternative='greater')
                                    print(f"{tissue} - {genetype} ({cistrans_counts_dict[genetype]}) - enrichment: {enrichment} - pval_binom {pval_binom}")
                                    enrichment_genetypes.append(CisTrans_Type(tissue=tissue, genetype=genetype, \
                                                                            ncistrans=TOTAL_CIS_TRANS, \
                                                                            ntype=cistrans_counts_dict[genetype], \
                                                                            frac_cis=frac_cis, \
                                                                            frac_cistrans=frac_cistrans, \
                                                                            enrichment=enrichment, \
                                                                            binom=pval_binom))

            targets_enrichment_outfile = os.path.join(baseoutdir,f"ciseqtl_target_enrichment_"+out_suffix)
            # if use_LD:
            #     targets_enrichment_outfile = os.path.join(baseoutdir,f"ciseqtl_target_enrichment_pc_lncRNA_{gamma_suffix}_knn30_cut{pcutoff}_ldpruned.txt")
            # else:
            #     targets_enrichment_outfile = os.path.join(baseoutdir,f"ciseqtl_target_enrichment_pc_lncRNA_{gamma_suffix}_knn30_cut{pcutoff}.txt")
            with open(targets_enrichment_outfile, 'w') as outstream:
                outstream.write("tissue\tgenetype\tncistrans\tntype\tfrac_cis\tfrac_cistrans\tenrichment\tpval_binom\n")
                for e in enrichment_genetypes:
                    outstream.write(f"{e.tissue}\t{e.genetype}\t{e.ncistrans}\t{e.ntype}\t{e.frac_cis}\t{e.frac_cistrans}\t{e.enrichment}\t{e.binom}\n")