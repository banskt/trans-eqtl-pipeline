import numpy as np
import os
import sys
import re
import collections
from shutil import copyfile
import argparse

def myreplace(s):
    todelete = ["(", ")", "-"]
    for ch in todelete:
        s = s.replace(ch, "", 1)
    return s.replace("  ", " ")

def read_tissues(infile, plain=False):
    tissues = []
    descriptions = []
    with open(infile) as instream:
        for l in instream:
            if re.search("^#", l) or re.search("^\s+$", l):
                continue
            tissues.append(l.split("\t")[1].rstrip())
            descriptions.append(l.split("\t")[0].rstrip())
    if not plain:
        descriptions = [myreplace(d) for d in descriptions]
    return tissues, descriptions

def parse_args():

    parser = argparse.ArgumentParser(description='Run DHS enrichments')

    parser.add_argument('--indir',
                        type=str,
                        dest='indir',
                        metavar='DIR',
                        help='Input directory of TEJAAS results')

    parser.add_argument('--tissue-file',
                        type=str,
                        dest='tissue_file',
                        help='Tissue file')

    parser.add_argument('--datatype',
                        type=str,
                        dest='datatype',
                        help='gtex or fhs')

    parser.add_argument('--preprocs',
                        nargs='+',
                        dest='preprocs',
                        help='Expression type [raw|lasso|others..]')

    opts = parser.parse_args()
    return opts


opts = parse_args()
sourcedir = opts.indir #"/cbscratch/franco/trans-eqtl/dev-pipeline/gtex_v8_a2s_lncRNA/"
pcutoff = 5e-08
if opts.datatype == "fhs":
    for preproc in opts.preprocs:
        print(sourcedir, preproc)
        #preproc = "permnull_sb{:s}_knn{:s}".format(sb2, K)
        basedir = sourcedir+"/{:s}/tejaas/{:s}".format(opts.datatype, preproc)
        destdir = sourcedir+"/summary_{:g}/{:s}/tejaas/{:s}".format(pcutoff, opts.datatype, preproc)
        if os.path.exists(basedir):
            filelist = [x for x in os.listdir(basedir) if x.startswith("t")]
            for file in filelist:
                if os.path.exists(os.path.join(basedir, file)):
                    if not os.path.exists(destdir): os.makedirs(destdir)
                if file == "target_genes_{:g}.txt".format(pcutoff): 
                    destfile = "target_genes.txt"
                    copyfile(os.path.join(basedir, file), os.path.join(destdir, destfile))
                if file == "target_genes_knn_{:g}.txt".format(pcutoff): 
                    destfile = "target_genes_knn.txt"
                    copyfile(os.path.join(basedir, file), os.path.join(destdir, destfile))
                if file == "trans_eqtls_{:g}.txt".format(pcutoff): 
                    destfile = "trans_eqtls.txt"
                    copyfile(os.path.join(basedir, file), os.path.join(destdir, destfile))
            if os.path.exists(destdir):
                file = "snps_list.txt"
                copyfile(os.path.join(basedir, file), os.path.join(destdir, file))
                
if opts.datatype == "gtex_v8":
    tissuenames, descriptions = read_tissues(opts.tissue_file)
    for t in tissuenames:
        for preproc in opts.preprocs:
            print(sourcedir, t, preproc)
            #preproc = "permnull_sb{:s}_knn{:s}".format(sb2, K)
            basedir = sourcedir+"/gtex_v8-{:s}/tejaas/{:s}".format(t, preproc)
            destdir = sourcedir+"/summary_{:g}/{:s}/tejaas/{:s}".format(pcutoff, t, preproc)
            if os.path.exists(basedir):
                filelist = [x for x in os.listdir(basedir) if x.startswith("t")]
                for file in filelist:
                    if os.path.exists(os.path.join(basedir, file)):
                        if not os.path.exists(destdir): os.makedirs(destdir)
                    if file == "target_genes_{:g}.txt".format(pcutoff): 
                        destfile = "target_genes.txt"
                        copyfile(os.path.join(basedir, file), os.path.join(destdir, destfile))
                    if file == "target_genes_knn_{:g}.txt".format(pcutoff): 
                        destfile = "target_genes_knn.txt"
                        copyfile(os.path.join(basedir, file), os.path.join(destdir, destfile))
                    if file == "trans_eqtls_{:g}.txt".format(pcutoff): 
                        destfile = "trans_eqtls.txt"
                        copyfile(os.path.join(basedir, file), os.path.join(destdir, destfile))
                if os.path.exists(destdir):
                    file = "snps_list.txt"
                    copyfile(os.path.join(basedir, file), os.path.join(destdir, file))