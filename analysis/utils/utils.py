import copy
import numpy as np
import re
import pandas as pd
from functools import wraps
import time

def timeit(f):
    @wraps(f)
    def wrap(*args, **kw):
        ts = time.time()
        result = f(*args, **kw)
        te = time.time()
        print('{:s} took: {:.6f} seconds'.format(f.__name__, te-ts))
        return result
    return wrap

def myreplace(s):
    todelete = ["(", ")", "-"]
    for ch in todelete:
        s = s.replace(ch, "")
    return s.replace("  ", " ")

def read_tissues(infile):
    tissues = []
    descriptions = []
    with open(infile) as instream:
        for l in instream:
            if re.search("^#", l):
                continue
            tissues.append(l.split("\t")[1].rstrip())
            descriptions.append(l.split("\t")[0].rstrip())
    descriptions = [myreplace(d) for d in descriptions]
    return tissues, descriptions

def read_rocfile(infile):
    df = pd.read_table(infile, header=0)
    nsel = np.array(df.nsel.tolist())
    tpr = np.array(df.tpr.tolist())
    ppv = np.array(df.ppv.tolist())
    valids = np.array(df.valids.tolist())
    return nsel, tpr, ppv, valids

@timeit
def get_compatible_snp_dicts(dict1, dict2):
    k1  = set(dict1.keys())
    k2  = set(dict2.keys())
    intersection = k1 & k2

    ndict1 = dict()
    ndict2 = dict()
    for k in intersection:
        ndict1[k] = dict1[k]
        ndict2[k] = dict2[k]
    return ndict1, ndict2

# def get_compatible_snp_dicts_old(dict1, dict2):
#     k1  = list(dict1.keys())
#     k2  = list(dict2.keys())

#     ndict1 = copy.deepcopy(dict1)  # takes ~ 1.21 s
#     ndict2 = copy.deepcopy(dict2)

#     # see if snps in dict1 are in dict2
#     for k in k1:
#         val2 = ndict2.get(k, None)
#         if val2 == None:
#             del ndict1[k]

#     for k in k2:
#         val1 = ndict1.get(k, None)
#         if val1 == None:
#             del ndict2[k]

#     return ndict1, ndict2
