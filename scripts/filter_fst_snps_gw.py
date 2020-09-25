import sys, os
import numpy as np
import collections
import argparse
import re

def parse_args():

    parser = argparse.ArgumentParser(description='Filter SNPs by Fst.')

    parser.add_argument('--tissues',
                        type=str,
                        dest='tissue_file',
                        metavar='FILE',
                        help='tissue list file')

    parser.add_argument('--fst',
                        type=str,
                        dest='fst_file',
                        help='file with Fst values')

    parser.add_argument('--outdir',
                        type=str,
                        dest='outdir',
                        help='base output directory')
    
    parser.add_argument('--basedir',
                        type=str,
                        dest='basedir',
                        help='base inputput directory')


    opts = parser.parse_args()
    return opts

def allreplace(s):
    todelete = ["(", ")"]
    for ch in todelete:
        s = s.replace(ch, "")
    s = s.replace(" - ", " ")
    s = s.replace("  ", " ")
    s = s.replace(" ", "_")
    return s

def read_tissues_str(infile):
    tissues = []
    descriptions = []
    tstrings = []
    with open(infile) as instream:
        for l in instream:
            lsp = l.split("\t")
            if re.search("^#", l):
                continue
            tissues.append(lsp[1].rstrip())
            descriptions.append(lsp[0].rstrip())
            tstrings.append(lsp[2].rstrip())
    #tstrings = [partreplace(d) for d in descriptions]
    descriptions = [allreplace(d) for d in descriptions]
    return tissues, descriptions, tstrings

def load_fst(fst_file):
    print("Loading fst values")
    fst_dict = dict()
    for chrm in range(1, 23):
        fst_dict[chrm] = collections.defaultdict(lambda: False)
    with open(fst_file) as instream:    
        next(instream)
        for line in instream:
            arr = line.strip().split()
            if len(arr) > 0:
                chrom = int(arr[0])
                pos   = int(arr[1])
                fstval= float(arr[2])
                fst_dict[chrom][pos] = fstval
    return fst_dict
if __name__ == '__main__':

    opts = parse_args()

    chroms = [str(x) for x in range(1,23)]
    fst_dict = load_fst(opts.fst_file)
    tissues, descriptions, tstrings = read_tissues_str(opts.tissue_file)

    if not os.path.exists(opts.outdir): os.makedirs(opts.outdir)

    genefiles = ["target_genes.txt", "target_genes_knn.txt"]

    for tissue in tissues:
        not_fst = 0
        reject_fst = 0
        teqtls_list = collections.defaultdict(lambda: False)
        tissue_outdir = os.path.join(opts.outdir, tissue)

        # Filter trans-eQTL list
        trans_file    = os.path.join(opts.basedir, tissue, "trans_eqtls.txt")
        if os.path.exists(trans_file):
            if not os.path.exists(tissue_outdir): os.makedirs(tissue_outdir)

            with open(os.path.join(tissue_outdir, "trans_eqtls.txt"), 'w') as outstream:
                with open(trans_file) as instream:
                    header = instream.readline()
                    outstream.write(header)
                    for line in instream:
                        arr = line.strip().split()
                        varid = arr[0]
                        chrom = int(varid.split("_")[0][3:])
                        pos   = int(varid.split("_")[1])
                        if pos in fst_dict[chrom]:
                            if fst_dict[chrom][pos] >= 0.3:
                                outstream.write(line)
                            else:
                                reject_fst += 1
                        else:
                            not_fst += 1
                            # print(f"{varid} not in Fst list")
            print(f"{tissue}: {not_fst} trans-eqtls not found in Fst list")
            print(f"{tissue}: {reject_fst} trans-eqtls filtered out by Fst")

            # Filter target gene list
            for genefile in genefiles:
                targets = os.path.join(opts.basedir, tissue, genefile)
                with open(targets) as instream:
                    with open(os.path.join(tissue_outdir, genefile), 'w') as outstream:
                        header = instream.readline()
                        outstream.write(header)
                        for line in instream:
                            varid = line.strip().split()[0]
                            if teqtls_list[varid]:
                                outstream.write(line)
        else:
            print(f"{tissue}: results file not found")