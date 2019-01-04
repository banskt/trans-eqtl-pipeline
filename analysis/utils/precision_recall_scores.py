import numpy as np

def confusion_matrix(result):
    keepres = [x for x in result if not np.isnan(x.stat)]
    nitems = len(keepres)
    ypred = np.array([x.stat for x in keepres])
    ytrue = np.array([x.causality for x in keepres])

    pos = len([x for x in result if x.causality == 1])
    neg = len(result) - pos

    tp = 0	# True positives
    fp = 0	# False positives

    tplist = list()
    fplist = list()

    if nitems > 0:
        alpha = np.max(ypred) + 1.0 # set threshold above the maximum value of predictions
        isort = np.argsort(ypred)[::-1]
    
    for j in range(nitems):
        if not ypred[isort[j]] == alpha:
            tplist.append(tp)
            fplist.append(fp)
            alpha = ypred[isort[j]]
        if ytrue[isort[j]] == 1:
            tp += 1
        else:
            fp += 1

    fplist.append(fp)
    tplist.append(tp)

    #tpr    = np.array([x / pos if pos > 0 else 0 for x in tplist]) 				# TPR = Recall = TP / Positives
    #fpr    = np.array([x / neg if neg > 0 else 0 for x in fplist]) 				# FPR = FP / Negatives
    #ppv    = np.array([x[0] / sum(x) if sum(x) > 0 else 1 for x in zip(tplist, fplist)]) 	# PPV = Precision = TP / (TP + FP)
    #fdr    = np.array([x[1] / sum(x) if sum(x) > 0 else 0 for x in zip(tplist, fplist)]) 	# FDR = FP / (TP + FP)
    tparr  = np.array(tplist)                                                                   # Number of true positives
    nsel   = np.array([sum(x) for x in zip(tplist, fplist)]) 					# Number of y selected at each threshold

    #return fpr, tpr, ppv, nsel, fdr
    return tparr, nsel
