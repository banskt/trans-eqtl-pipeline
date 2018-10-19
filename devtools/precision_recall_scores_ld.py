import numpy as np

def confusion_matrix(ldresult, ldcut = 1.0):
    keepres = [x for x in ldresult if not np.isnan(x.stat)]
    nitems = len(keepres)
    ypred = np.array([x.stat for x in keepres])
    ytrue = np.array([x.causality for x in keepres])

    pos = len([x for x in ldresult if x.causality == 1])
    neg = len(ldresult) - pos

    tp = 0	# True positives
    fp = 0	# False positives
    tpld = 0	# True positives (if in strong LD)
    fpld = 0    # False positives (not in LD with actual TP)

    tplist = list()
    fplist = list()
    tpldlist = list()
    fpldlist = list()

    if nitems > 0:
        alpha = np.max(ypred) + 1.0 # set threshold above the maximum value of predictions
        isort = np.argsort(ypred)[::-1]
    
    for j in range(nitems):
        if not ypred[isort[j]] == alpha:
            tplist.append(tp)
            fplist.append(fp)
            tpldlist.append(tpld)
            fpldlist.append(fpld)
            alpha = ypred[isort[j]]
        if ytrue[isort[j]] == 1:
            tp += 1
        else:
            fp += 1
        if ldresult[isort[j]].ld > ldcut:
            tpld += 1
        else:
            fpld += 1

    fplist.append(fp)
    tplist.append(tp)
    tpldlist.append(tpld)
    fpldlist.append(fpld)

    tpr    = np.array([x / pos if pos > 0 else 0 for x in tplist]) 				# TPR = Recall = TP / Positives
    ldtpr  = np.array([x / pos if pos > 0 else 0 for x in tpldlist])
    fpr    = np.array([x / neg if neg > 0 else 0 for x in fplist]) 				# FPR = FP / Negatives
    ppv    = np.array([x[0] / sum(x) if sum(x) > 0 else 1 for x in zip(tplist, fplist)]) 	# PPV = Precision = TP / (TP + FP)
    ldppv  = np.array([x[0] / sum(x) if sum(x) > 0 else 1 for x in zip(tpldlist, fpldlist)])
    fdr    = np.array([x[1] / sum(x) if sum(x) > 0 else 0 for x in zip(tplist, fplist)]) 	# FDR = FP / (TP + FP)
    nsel   = np.array([sum(x) for x in zip(tplist, fplist)]) 					# Number of y selected at each threshold
    ldnsel = np.array([sum(x) for x in zip(tpldlist, fpldlist)])

    return fpr, tpr, ppv, nsel, fdr, ldtpr, ldppv, ldnsel
