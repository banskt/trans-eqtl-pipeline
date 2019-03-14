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
# mpl_stylesheet.banskt_presentation(fontfamily = 'system')

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


def validate(testdict, valdict, smax = -np.log10(0.05), nmax = None, empirical = False):
    if empirical:
        if nmax is None:
            totsnps = len(list(testdict.keys()))
            nmax = int(np.round(totsnps/20))
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
    datalist = list()
    for key, value in data.items():
        datalist.append(value)
    return datalist

chrms = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
# datasets = ['gtex-hlv', 'gtex-ms']
# datasets = ['gtex-wb', 'gtex-sse', 'gtex-thy', 'gtex-bca', 'gtex-bce', 'gtex-aa']
# datasets = ["gtex-ban", "gtex-bceh"] #, "gtex-bco", "gtex-bfr", "gtex-bhi", "gtex-bhy","gtex-bnu", "gtex-bpu", "gtex-bca"]
datasets = ["gtex-ban", "gtex-bceh", "gtex-bco", "gtex-bfr", "gtex-bhi", "gtex-bhy", \
"gtex-bnu", "gtex-bpu", "gtex-bca", "gtex-bce", "gtex-as", \
"gtex-av", "gtex-ag", "gtex-ac", "gtex-at", "gtex-br", "gtex-ebv", \
"gtex-fib", "gtex-cols", "gtex-colt", "gtex-esog", "gtex-esom", "gtex-esomu", \
"gtex-haa", "gtex-liv", "gtex-lu", "gtex-ov", "gtex-pan", "gtex-pit", "gtex-pro", \
"gtex-snse", "gtex-si", "gtex-spl", "gtex-sto", "gtex-ut", "gtex-va", \
"gtex-aa", "gtex-hlv", "gtex-ms", "gtex-sse", "gtex-tes", "gtex-thy", "gtex-wb"]

# methods = ['tejaas_maf', 'tejaas_perm', 'cpma', 'matrixeqtl', 'matrixeqtl_fdr',
#            'tejaas_rand_maf', 'tejaas_rand_perm', 'cpma_rand', 'matrixeqtl_rand']

methods = ['tejaas_perm', 'matrixeqtl', 'matrixeqtl_fdr'] #, 'tejaas_rand_perm', 'matrixeqtl_fdr', 'matrixeqtl_rand']

input_dir = '/cbscratch/franco/trans-eqtl/dev-pipeline/lmcorrected'
combi = [x for x in itertools.combinations(datasets,2)]

# Load all datasets first
data = collections.defaultdict(dict)
for dataset in datasets:
    for key in methods:
        if input_dir.endswith('uncorrected'):
            data[dataset][key] = get_dict_for_no_crxn(key, input_dir, chrms, dataset)
        else:
            data[dataset][key] = get_dict_for_method(key, input_dir, chrms, dataset)
        if key.endswith("_sp"):
            print("going with ", dataset, key)
            data[dataset][key] = update_sparse_dict(data[dataset][key.rstrip("_sp")].copy(), data[dataset][key])

# generate all pairs of test x crossvalidations (maybe to much memory?)
res = collections.defaultdict(dict)
tissue_pairs = []
for combination in combi:
    tissue_pair = combination[0]+"_"+combination[1]
    tissue_pairs.append(tissue_pair)
    for key in methods:
        print(key, tissue_pair)
        testset     = data[combination[0]][key]
        crossvalset = data[combination[1]][key]
        roc_file = tissue_pair+"."+key+".roc.txt"
        if not os.path.exists(roc_file):
            testdict, valdict = utils.get_compatible_snp_dicts(testset, crossvalset)
            valresult = validate(testdict, valdict, empirical=True)
            nsel, tpr, ppv, valids = roc.confusion_matrix(valresult)
            with open("newrun/"+roc_file, 'w') as outstream2:
                outstream2.write("nsel\ttpr\tppv\tvalids\n")
                for i in range(len(nsel)):
                    outstream2.write("{:g}\t{:g}\t{:g}\t{:g}\n".format(nsel[i], tpr[i], ppv[i], valids[i]))
        else:
            print("File exists:",roc_file)
            nsel, tpr, ppv, valids = utils.read_rocfile(roc_file)
        scaled_nsel = nsel / max(nsel)
        auc = metrics.auc(scaled_nsel, tpr)
        with open("newrun/"+key+".auc.txt", 'a') as outstream:
            outstream.write("{:s}\t{:s}\t{:g}\n".format(tissue_pair.split("_")[0], tissue_pair.split("_")[1], auc))