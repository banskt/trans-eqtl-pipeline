import copy
import numpy as np
import re
import pandas as pd

def read_tissues(infile):
    tissues = []
    descriptions = []
    with open(infile) as instream:
        for l in instream:
            if re.search("^#", l):
                continue
            tissues.append(l.split("\t")[1].rstrip())
            descriptions.append(l.split("\t")[0].rstrip().replace(" ", "_"))
    return tissues, descriptions

def read_rocfile(infile):
    df = pd.read_table(infile, header=0)
    nsel = np.array(df.nsel.tolist())
    tpr = np.array(df.tpr.tolist())
    ppv = np.array(df.ppv.tolist())
    valids = np.array(df.valids.tolist())
    return nsel, tpr, ppv, valids

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
