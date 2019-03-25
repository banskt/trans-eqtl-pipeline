import numpy as np
import collections
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import os
import sys
import itertools 
sys.path.append('../')
from utils import load_results
from utils import precision_recall_scores as roc
from utils import mpl_stylesheet
from utils import utils
from sklearn import metrics
import time
# mpl_stylesheet.banskt_presentation(fontfamily = 'system')
import argparse

def parse_args():

    parser = argparse.ArgumentParser(description='Convert rpkm for a given tissue into inverse quantile normalized gene expression')

    parser.add_argument('--method',
                        type=str,
                        dest='method',
                        help='options: tejaas_perm, matrixeqtl, matrixeqtl_fdr')

    parser.add_argument('--sb',
                        type=str,
                        dest='sigmabeta',
                        metavar='SIGMA',
                        help='sigma beta for tejaas')

    opts = parser.parse_args()
    return opts



def get_dict_for_method(method, input_dir, chrms, tissue, sb = '0.01'):
    chrmdicts = [dict() for c in chrms]
    datadir = os.path.join(input_dir, tissue)
    print("Reading {:s} for {:s}".format(tissue, method))
    for i, chrm in enumerate(chrms):
        if method == 'tejaas_maf':
            sb = '0.001'
            filepath = os.path.join(datadir, 'tejaas', 'mafnull_sb{:s}'.format(sb), 'chr{:d}'.format(chrm), 'rr.txt')
            chrmdicts[i] = load_results.tejaas(filepath)
        elif method == 'tejaas_rand_maf':
            sb = '0.001'
            filepath = os.path.join(datadir, 'tejaas_rand', 'mafnull_sb{:s}'.format(sb), 'chr{:d}'.format(chrm), 'rr.txt')
            chrmdicts[i] = load_results.tejaas(filepath)
        elif method == 'tejaas_perm':
            filepath = os.path.join(datadir, 'tejaas', 'permnull_sb{:s}'.format(sb), 'chr{:d}'.format(chrm), 'rr.txt')
            print("Loading ", filepath)
            chrmdicts[i] = load_results.tejaas(filepath)
        elif method == 'tejaas_perm_sp':
            filepath = os.path.join(datadir, 'tejaas', 'permnull_sb{:s}'.format(sb), 'chr{:d}'.format(chrm), 'rr_it1.txt')
            print("Loading ", filepath)
            chrmdicts[i] = load_results.tejaas(filepath)
        elif method == 'tejaas_rand_perm':
            filepath = os.path.join(datadir, 'tejaas_rand', 'permnull_sb{:s}'.format(sb), 'chr{:d}'.format(chrm), 'rr.txt')
            print("Loading ", filepath)
            chrmdicts[i] = load_results.tejaas(filepath)
        elif method == 'tejaas_rand_perm_sp':
            filepath = os.path.join(datadir, 'tejaas_rand', 'permnull_sb{:s}'.format(sb), 'chr{:d}'.format(chrm), 'rr_it1.txt')
            print("Loading ", filepath)
            chrmdicts[i] = load_results.tejaas(filepath)
        elif method == 'cpma':
            filepath = os.path.join(datadir, 'tejaas', 'jpa', 'chr{:d}'.format(chrm), 'jpa.txt')
            chrmdicts[i] = load_results.jpa(filepath)
        elif method == 'cpma_rand':
            filepath = os.path.join(datadir, 'tejaas_rand', 'jpa', 'chr{:d}'.format(chrm), 'jpa.txt')
            chrmdicts[i] = load_results.jpa(filepath)
        elif method == 'matrixeqtl':
            filepath = os.path.join(datadir, 'matrixeqtl', 'chr{:d}'.format(chrm), 'trans_eqtl.txt')
            print("Loading ", filepath)
            chrmdicts[i] = load_results.matrixeqtl(filepath)
        elif method == 'matrixeqtl_fdr':
            filepath = os.path.join(datadir, 'matrixeqtl', 'chr{:d}'.format(chrm), 'trans_eqtl.txt')
            chrmdicts[i] = load_results.matrixeqtl_fdr(filepath)
        elif method == 'matrixeqtl_rand':
            filepath = os.path.join(datadir, 'matrixeqtl_rand', 'chr{:d}'.format(chrm), 'trans_eqtl.txt')
            chrmdicts[i] = load_results.matrixeqtl(filepath)

    res = dict()
    for d in chrmdicts:
        for k, v in d.items():
            res[k] = v
    return res

def update_sparse_dict(old_dict, new_dict):
    for key in new_dict.keys():
        if old_dict.get(key, False):
            old_dict[key] = new_dict[key]
        else:
            print("SNP {:s} not present".format(key))
    return old_dict

def get_dict_for_no_crxn(method, input_dir, chrms, tissue, sb = '0.01'):
    chrmdicts = [dict() for c in chrms]
    datadir = os.path.join(input_dir, tissue)
    print("Reading {:s} for {:s}".format(tissue, method))
    for i, chrm in enumerate(chrms):
        if method == 'matrixeqtl':
            filepath = os.path.join(datadir, 'gtex_MatrixEQTL_chr{:d}.transout'.format(chrm))
            chrmdicts[i] = load_results.matrixeqtl(filepath)
        elif method == 'matrixeqtl_fdr':
            filepath = os.path.join(datadir, 'gtex_MatrixEQTL_chr{:d}.transout'.format(chrm))
            chrmdicts[i] = load_results.matrixeqtl_fdr(filepath)

    res = dict()
    for d in chrmdicts:
        for k, v in d.items():
            res[k] = v
    return res

INFO_FIELDS = ['rsid', 'stat', 'causality']
class ValidateResult(collections.namedtuple('_ValidateResult', INFO_FIELDS)):
    __slots__ = ()

def validate(testdict, valdict, smax = -np.log10(0.05), nmax = None, empirical = False, empirical_percent = 5):
    if empirical:
        percent_fraction = 100/empirical_percent
        if nmax is None:
            totsnps = len(list(testdict.keys()))
            nmax = int(np.round(totsnps/percent_fraction))
            print("totsnps:", totsnps, " nmax:", nmax)
        sorted_snps = [x[0] for x in sorted(valdict.items(), key=lambda x: x[1], reverse = True)]
        topsnps = sorted_snps[:nmax]
    else:
        topsnps = [k for k,v in valdict.items() if v > smax]
    data = dict()
    for key, value in testdict.items():
        data[key] = ValidateResult(rsid = key, stat = value, causality = 0)
    for key in topsnps:
        if key in data:
            data[key] = ValidateResult(rsid = key, stat = data[key].stat, causality = 1)
    datalist = list(data.values())
    return datalist

def calc_auc(nsel, tpr, n=1000):
    scaled_nsel = nsel / max(nsel)
    auc = metrics.auc(scaled_nsel, tpr)
    ix = nsel < n
    xscaled  = np.append(scaled_nsel[ ix ], n/max(nsel))
    ytpr = np.append(tpr[ ix ], tpr[ np.sum(ix)])
    auc1k = metrics.auc(xscaled, ytpr)
    return auc, auc1k

@utils.timeit
def validate_and_roc(testset, crossvalset, tissue_pair, options, roc_file, auc_file):
    if not os.path.exists(roc_file):
        valresult = validate(testdict, valdict, empirical=options["empirical"], empirical_percent=options["empirical_percent"])
        nsel, tpr, ppv, valids = roc.confusion_matrix(valresult)
        with open(roc_file, 'w') as outstream2:
            outstream2.write("nsel\ttpr\tppv\tvalids\n")
            for i in range(len(nsel)):
                outstream2.write("{:g}\t{:g}\t{:g}\t{:g}\n".format(nsel[i], tpr[i], ppv[i], valids[i]))
    else:
        print("File exists:",roc_file)
        nsel, tpr, ppv, valids = utils.read_rocfile(roc_file)
    auc, auc1k = calc_auc(nsel, tpr, n=1000)
    with open(auc_file, 'a') as outstream:
        outstream.write("{:s}\t{:s}\t{:g}\t{:g}\n".format(tissue_pair.split("_")[0], tissue_pair.split("_")[1], auc, auc1k))

if __name__ == '__main__':
    opts = parse_args()

    if opts.method not in ['tejaas_maf', 'tejaas_perm', 'cpma', 'matrixeqtl', 'matrixeqtl_fdr',
                            'tejaas_rand_maf', 'tejaas_rand_perm', 'cpma_rand', 'matrixeqtl_rand']:
        print("ERROR METHOD", opts.method)
        raise

    ts = time.time()
    chrms = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
    # sbs = ["0.01", "0.05"]  # sigma_betas
    sbs = [opts.sigmabeta]
    eps = [5, 1, 0.05]        # empirical_percents

    tissue_file = "/usr/users/fsimone/trans-eqtl-pipeline/analysis/plots/tissues.txt"
    tissues, descriptions = utils.read_tissues(tissue_file)
    datasets = ["gtex-"+t for t in tissues]

    # methods = ['tejaas_maf', 'tejaas_perm', 'cpma', 'matrixeqtl', 'matrixeqtl_fdr',
    #            'tejaas_rand_maf', 'tejaas_rand_perm', 'cpma_rand', 'matrixeqtl_rand']

    methods = [opts.method]
    
    input_dir = '/cbscratch/franco/trans-eqtl/dev-pipeline/lmcorrected'
    outdirbase = "/cbscratch/franco/trans-eqtl/analysis/data/"
    combi = [x for x in itertools.combinations(datasets,2)]

    # Load all datasets first
    allmethods = list()
    data = collections.defaultdict(dict)
    for dataset in datasets:
        for key in methods:
            if key.startswith("tejaas"):
                for sb in sbs:
                    mkey = key + "_" + sb
                    print("Loading",mkey)
                    if input_dir.endswith('uncorrected'):
                        data[dataset][mkey] = get_dict_for_no_crxn(key, input_dir, chrms, dataset, sb)
                    else:
                        data[dataset][mkey] = get_dict_for_method(key, input_dir, chrms, dataset, sb)
                    if mkey not in allmethods: allmethods.append(mkey) 
                    # broken...
                    # if key.endswith("_sp"):
                    #     print("going with ", dataset, mkey)
                    #     data[dataset][key] = update_sparse_dict(data[dataset][key.rstrip("_sp")].copy(), data[dataset][key])
            else:
                print("Loading",key)
                if input_dir.endswith('uncorrected'):
                    data[dataset][key] = get_dict_for_no_crxn(key, input_dir, chrms, dataset)
                else:
                    data[dataset][key] = get_dict_for_method(key, input_dir, chrms, dataset)
                if key not in allmethods: allmethods.append(key) 

    #init objects
    # pool = mp.Pool(16)

    res = collections.defaultdict(dict)
    tissue_pairs = []
    for combination in combi:
        jobs_val = []
        tissue_pair = combination[0]+"_"+combination[1]
        tissue_pairs.append(tissue_pair)
        for key in allmethods:
            testset     = data[combination[0]][key]
            crossvalset = data[combination[1]][key]
            print(combi, key)
            testdict, valdict = utils.get_compatible_snp_dicts(testset, crossvalset)
            for ep in eps:
                print(key, ep, tissue_pair)
                outdir = os.path.join(outdirbase,key,str(ep))
                os.makedirs(outdir, exist_ok=True)
                options = dict()
                options["smax"] = -np.log10(0.05)
                options["empirical"] = True
                options["empirical_percent"] = ep
                roc_file = os.path.join(outdir, tissue_pair+"."+key+".roc.txt")
                auc_file = os.path.join(outdir, key+".auc.txt")
                validate_and_roc(testdict, valdict, tissue_pair, options, roc_file, auc_file)
                #jobs_val.append( pool.apply_async( validate_and_roc, (testdict, valdict, tissue_pair, options, roc_file, auc_file)))

    # for job in jobs_val:
    #     job.get()

    # pool.close()
    te = time.time()
    print('Validation took: {:.6f} seconds'.format(te-ts))