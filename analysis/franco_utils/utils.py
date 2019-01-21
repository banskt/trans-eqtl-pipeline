import copy
import mpmath
import numpy as np
import os 
import random
import operator
import collections

# OPTIONS_FIELDS = ['tejaas_method', 'pval_thres', 'zoom', 'zoom_percent']
# class Options(collections.namedtuple('_Options', OPTIONS_FIELDS)):
#     __slots__ = ()

# decimal places
mpmath.mp.dps = 500
def pval(x): return mpmath.log10(1 - 0.5 * (1 + mpmath.erf(x/mpmath.sqrt(2))))

def load_tejaas_jpa_results(chr_list, input_dir):
    tejaas_dict = dict()
    for chrm in chr_list:
        print( "Reading TEJAAS chr{:d} from {:s}".format(chrm, input_dir))
        dirc = os.path.join(input_dir, "chr{:d}".format(chrm))
        pths = [os.path.join(dirc,path) for path in os.listdir(dirc) if "_jpa.txt" in path]
        for pth in pths:
            pvals= list()
            l = open(str(pth),"r").readlines()
            for line in l[1:]:    
                arr      = line.strip().split("\t")
                rsid     = arr[0]
                jpascore = float(arr[1])
                tejaas_dict[rsid] = jpascore
    return tejaas_dict

def load_tejaas_results(chr_list, input_dir):
    tejaas_dict = dict()
    for chrm in chr_list:
        print( "Reading TEJAAS chr{:d} from {:s}".format(chrm, input_dir))
        dirc = os.path.join(input_dir, "chr{:d}".format(chrm))
        pths = [os.path.join(dirc,path) for path in os.listdir(dirc) if "_rr.txt" in path]
        for pth in pths:
            pvals= list()
            l = open(str(pth),"r").readlines()
            for line in l[1:]:    
                arr   = line.strip().split("\t")
                rsid  = arr[0]
                P     = float(arr[5])
                Q     = float(arr[2])
                Mu    = float(arr[3])
                Sigma = float(arr[4])
                pvalue= np.log10(P) if P!=0 else pval((Q-Mu)/Sigma)
                tejaas_dict[rsid] = pvalue 
    return tejaas_dict

def load_matrixeqtl_results(chr_list, input_file):
    matrix_dict = {}
    for chrm in chr_list:
        print( "Reading matrixEQTL chr{:d} from {:s}".format(chrm, input_file))
        results_file = os.path.join(input_file.format(chrm))
        with open(results_file) as instream:
            _ = instream.readline()
            for line in instream:
                arr  = line.rstrip().split("\t")
                rsid = arr[0]
                FDR  = float(arr[5])
                if rsid not in matrix_dict:
                    matrix_dict[rsid] = np.log10(FDR)
    return matrix_dict

def get_compatible_snp_dicts(dict1, dict2):
    k1  = list(dict1.keys())
    k2  = list(dict2.keys())

    ndict1 = copy.deepcopy(dict1)  # takes ~ 1.21 s
    ndict2 = copy.deepcopy(dict2)    

    # see if snps in dict1 are in dict2
    for k in k1:
        val2 = ndict2.get(k, None) 
        if val2 == None:
            del ndict1[k]

    for k in k2:
        val1 = ndict1.get(k, None) 
        if val1 == None:
            del ndict2[k]

    return ndict1, ndict2


def get_sorted_dict_tuples(mydict):
    isReversed = False
    v_keys = list(mydict.keys())
    n_snps = len(v_keys)
    my_tuples = sorted(my_dict.items(), key=operator.itemgetter(1),reverse=isReversed)
    my_rsids = [my_tuples[i][0] for i in range(n_snps)]
    my_pvals = [my_tuples[i][1] for i in range(n_snps)]
    return my_rsids, my_pvals

def get_empirical_replication_sizes(test_dict, validation_dict, title, pval_thres=np.log10(0.05)):        
#     pval_thres = np.log10(0.05)
    
    rsids1  = list(test_dict.keys())
    rsids2  = list(validation_dict.keys())

    n_snps1 = len(rsids1)
    n_snps2 = len(rsids2)

#     common_snps = [snp for snp in rsids1 if validation_dict.get(snp, None) != None ]
    compat_test_dict, compat_validation_dict = get_compatible_snp_dicts(test_dict, validation_dict)
    compat_rsids = list(compat_test_dict.keys())
    common_n_snps = len(compat_rsids)

    validation_rsids, validation_pvals = get_sorted_dict_tuples(validation_dict)

    empirical005 = np.round(n_snps/20)
    label_dict = dict()
    ntrue = 0
    for i, k in enumerate(validation_rsids):
        ntrue += 1
        label_dict[k] = True if i < empirical005 else False

    print("ntrue: ", ntrue)


    test_tuples = sorted(test_dict.items(), key=operator.itemgetter(1),reverse=isReversed)
    test_rsids = [test_tuples[i][0] for i in range(n_snps)]
    test_pvals = [test_tuples[i][1] for i in range(n_snps)]





    signif_snps1 = [snp for snp in rsids1 if test_dict[snp] < pval_thres]
    signif_snps2 = [snp for snp in rsids2 if validation_dict[snp] < pval_thres]

    nsig1 = len(signif_snps1)
    nsig2 = len(signif_snps2)
    
    compat_signig_snps1 = [snp for snp in compat_rsids if compat_test_dict[snp] < pval_thres]
    compat_signig_snps2 = [snp for snp in compat_rsids if compat_validation_dict[snp] < pval_thres]

    common_signif_snps = [snp for snp in signif_snps1 if snp in signif_snps2]

    print("Algorithm      |Nº snps1|Nº snps2| common |  %   |Nº sig1 |Nº sig2 |NºC sig1|NºC sig2| common |   %  ")

    if n_snps1 > 0 and n_snps2 > 0:
        int_ratio = 100*common_n_snps/min(n_snps1, n_snps2)
    else:
        int_ratio = 0
    if nsig1 > 0 and nsig2 > 0:
        sig_ratio = 100*len(common_signif_snps)/min(nsig1, nsig2)
    else:
        sig_ratio = 0
    print("{:15s}|{:8d}|{:8d}|{:8d}|{:6.2f}|{:8d}|{:8d}|{:8d}|{:8d}|{:8d}|{:6.2f}".format(title, n_snps1,\
        n_snps2, common_n_snps, int_ratio, nsig1, nsig2, \
        len(compat_signig_snps1), len(compat_signig_snps2), \
        len(common_signif_snps), sig_ratio))


def get_replication_sizes(test_dict, validation_dict, title, pval_thres=np.log10(0.05)):        
#     pval_thres = np.log10(0.05)
    
    rsids1  = list(test_dict.keys())
    rsids2  = list(validation_dict.keys())

    n_snps1 = len(rsids1)
    n_snps2 = len(rsids2)

#     common_snps = [snp for snp in rsids1 if validation_dict.get(snp, None) != None ]
    compat_test_dict, compat_validation_dict = get_compatible_snp_dicts(test_dict, validation_dict)
    compat_rsids = list(compat_test_dict.keys())
    common_n_snps = len(compat_rsids)

    signif_snps1 = [snp for snp in rsids1 if test_dict[snp] < pval_thres]
    signif_snps2 = [snp for snp in rsids2 if validation_dict[snp] < pval_thres]

    nsig1 = len(signif_snps1)
    nsig2 = len(signif_snps2)
    
    compat_signig_snps1 = [snp for snp in compat_rsids if compat_test_dict[snp] < pval_thres]
    compat_signig_snps2 = [snp for snp in compat_rsids if compat_validation_dict[snp] < pval_thres]

    common_signif_snps = [snp for snp in signif_snps1 if snp in signif_snps2]

    print("Algorithm      |Nº snps1|Nº snps2| common |  %   |Nº sig1 |Nº sig2 |NºC sig1|NºC sig2| common |   %  ")

    if n_snps1 > 0 and n_snps2 > 0:
        int_ratio = 100*common_n_snps/min(n_snps1, n_snps2)
    else:
        int_ratio = 0
    if nsig1 > 0 and nsig2 > 0:
        sig_ratio = 100*len(common_signif_snps)/min(nsig1, nsig2)
    else:
        sig_ratio = 0
    print("{:15s}|{:8d}|{:8d}|{:8d}|{:6.2f}|{:8d}|{:8d}|{:8d}|{:8d}|{:8d}|{:6.2f}".format(title, n_snps1,\
        n_snps2, common_n_snps, int_ratio, nsig1, nsig2, \
        len(compat_signig_snps1), len(compat_signig_snps2), \
        len(common_signif_snps), sig_ratio))

def get_validation_curve_jpa(test_dict, validation_dict, thres=20, randomize=False):    
    n_snps = len(validation_dict.keys())
    isReversed = True # descending order
    validation_tuples = sorted(validation_dict.items(), key=operator.itemgetter(1),reverse=isReversed)
    validation_rsids = [validation_tuples[i][0] for i in range(n_snps)]
    validation_pvals = [validation_tuples[i][1] for i in range(n_snps)]

    # test_pvals = [test_dict[e] for e in test_dict.keys()]
    
    # print(np.min(validation_pvals), np.max(validation_pvals))
    # print(np.min(test_pvals), np.max(test_pvals))
    if randomize:
        random.shuffle(validation_rsids)
    
    toplot = []
    check_x = []
    positives = []
        
    i = 0 
    while(i < n_snps):
        try:
            if(test_dict[validation_rsids[i]] > thres):
                positives.append(i)
        except:
            pass
        i = i + 1
        while(i < n_snps and validation_pvals[i] == validation_pvals[i-1]):
            try:
                if(test_dict[validation_rsids[i]] > thres):
                    positives.append(i)
            except:
                pass
            i = i + 1
        check_x.append(i)
        toplot.append(len(positives))
    return check_x, toplot

def get_validation_curve(test_dict, validation_dict, pval_thres=0.05, islog=True, randomize=False):    
    n_snps = len(validation_dict.keys())
    isReversed = True
    if islog:
        isReversed = False
    validation_tuples = sorted(validation_dict.items(), key=operator.itemgetter(1),reverse=isReversed)
    validation_rsids = [validation_tuples[i][0] for i in range(n_snps)]
    validation_pvals = [validation_tuples[i][1] for i in range(n_snps)]

    # test_pvals = [test_dict[e] for e in test_dict.keys()]
    
    # print(np.min(validation_pvals), np.max(validation_pvals))
    # print(np.min(test_pvals), np.max(test_pvals))
    if randomize:
        random.shuffle(validation_rsids)
    
    toplot = []
    check_x = []
    positives = []
    if islog:
        pval_thres = np.log10(pval_thres)
        
    i = 0 
    while(i < n_snps):
        try:
            if(test_dict[validation_rsids[i]] < pval_thres):
                positives.append(i)
        except:
            pass
        i = i + 1
        while(i < n_snps and validation_pvals[i] == validation_pvals[i-1]):
            try:
                if(test_dict[validation_rsids[i]] < pval_thres):
                    positives.append(i)
            except:
                pass
            i = i + 1
        check_x.append(i)
        toplot.append(len(positives))
    return check_x, toplot

def get_roc_curve(statistic, label):
    print("getting ROC curve")
    x_vals = []
    i_vals = []
    n_vals = []
    recall = []
    ppv    = []
    total_TP = np.sum(label)
    TP = 0
    FP = 0
    n  = 0
    if label[0]:
        TP += 1
    else:
        FP += 1

    x_vals.append(statistic[0])
    i_vals.append(0)
    n_vals.append(n)
    recall.append(TP/total_TP)
    ppv.append(TP/(TP+FP))  
    for i in range(1, len(statistic)):
        if statistic[i] > statistic[i-1]:
            n += 1
            x_vals.append(statistic[i])
            n_vals.append(n)
            i_vals.append(i)
            recall.append(TP/total_TP)
            ppv.append(TP/(TP+FP))
        if label[i]:
            TP += 1
        else:
            FP += 1
    # x_vals.append(statistic[-1])
    # n_vals.append(n+1)
    # i_vals.append(len(statistic))
    # recall.append(TP/total_TP)
    # ppv.append(TP/(TP+FP))
    print("roc curve finished")

    return np.array(x_vals), np.array(i_vals), np.array(n_vals), np.array(recall), np.array(ppv)

def evaluate_replication_jpa(test_dict, validation_dict, randomize=False, thres = 20):
    isReversed = True  

    v_keys = list(validation_dict.keys())
    n_snps = len(v_keys)

    label_dict = dict()
    for k in v_keys:
        label_dict[k] = True if validation_dict[k] >= thres else False

    test_tuples = sorted(test_dict.items(), key=operator.itemgetter(1),reverse=isReversed)
    test_rsids = [test_tuples[i][0] for i in range(n_snps)]
    test_pvals = [test_tuples[i][1] for i in range(n_snps)]
    
    if randomize:
        random.shuffle(test_rsids)
        print("Randomized")
    
    # sort labels by pvalue in test set (rsids are already sorted by pval)
    label = [label_dict[snp] for snp in test_rsids]
    x_vals, i_vals, n_vals, recall, ppv = get_roc_curve(-np.array(test_pvals), label)
    
    return x_vals, i_vals, n_vals, recall, ppv

def evaluate_replication(test_dict, validation_dict, islog=True, randomize=False, pval_thres = 0.05):
    isReversed = True

    if islog:
        isReversed = False
        pval_thres = np.log10(pval_thres)
    

    v_keys = list(validation_dict.keys())
    n_snps = len(v_keys)

    label_dict = dict()
    for k in v_keys:
        label_dict[k] = True if validation_dict[k] <= pval_thres else False

    test_tuples = sorted(test_dict.items(), key=operator.itemgetter(1),reverse=isReversed)
    test_rsids = [test_tuples[i][0] for i in range(n_snps)]
    test_pvals = [test_tuples[i][1] for i in range(n_snps)]
    
    if randomize:
        random.shuffle(test_rsids)
        print("Randomized")
    
    # sort labels by pvalue in test set (rsids are already sorted by pval)
    label = [label_dict[snp] for snp in test_rsids]
    x_vals, i_vals, n_vals, recall, ppv = get_roc_curve(test_pvals, label)
    
    return x_vals, i_vals, n_vals, recall, ppv

def evaluate_replication_empirical(test_dict, validation_dict, islog=True, randomize=False, pval_thres = 0.05):
    isReversed = True

    if islog:
        isReversed = False

    v_keys = list(validation_dict.keys())
    n_snps = len(v_keys)

    validation_tuples = sorted(validation_dict.items(), key=operator.itemgetter(1),reverse=isReversed)
    validation_rsids = [validation_tuples[i][0] for i in range(n_snps)]
    validation_pvals = [validation_tuples[i][1] for i in range(n_snps)]

    empirical005 = np.round(n_snps/20)
    label_dict = dict()
    for i, k in enumerate(validation_rsids):
        label_dict[k] = True if i < empirical005 else False

    test_tuples = sorted(test_dict.items(), key=operator.itemgetter(1),reverse=isReversed)
    test_rsids = [test_tuples[i][0] for i in range(n_snps)]
    test_pvals = [test_tuples[i][1] for i in range(n_snps)]

    if randomize:
        random.shuffle(test_rsids)
        print("Randomized")

    # sort labels by pvalue in test set (rsids are already sorted by pval)
    label = [label_dict[snp] for snp in test_rsids]
    x_vals, i_vals, n_vals, recall, ppv = get_roc_curve(test_pvals, label)

    return x_vals, i_vals, n_vals, recall, ppv

def get_ranks(my_dict, descending=True, randomize=False):
    my_keys = list(my_dict.keys())
    n_snps = len(my_keys)

    tuples = sorted(my_dict.items(), key=operator.itemgetter(1),reverse=descending)
    rsids = [tuples[i][0] for i in range(n_snps)]
    pvals = [tuples[i][1] for i in range(n_snps)]

    if randomize:
        random.shuffle(rsids)
        print("Randomized")

    x_rank_dict = dict()
    for i, k in enumerate(rsids):
        x_rank_dict[k] = i
    return x_rank_dict