#!/usr/bin/env python

''' Unfortunately there is no standard for gtf files.
    So, every version gets different function
'''

import os
import numpy as np
import gzip
from containers import GeneInfo
import re

def gencode_v12(filepath, feature = 'gene', trim=False, biotype=['protein_coding', 'lncRNA'], include_chrom = 0, include_chroms=['{:d}'.format(x + 1) for x in range(22)]):
    annotfile = os.path.realpath(filepath)
    geneinfo = list()
    lncRNA_list = ["macro_lncRNA","non_coding","bidirectional_promoter_lncRNA","3prime_overlapping_ncRNA","sense_overlapping","processed_transcript","sense_intronic","antisense","lincRNA", "miRNA"]
    mode = "v19"
    if re.search("v26", filepath):
        mode = "v26"

    if "lncRNA" in biotype:
        mybiotype = biotype + lncRNA_list
    else:
        mybiotype = biotype
    if annotfile.endswith("affy"):
        geneinfo = affy_exon_chip(annotfile, include_chrom = include_chrom, include_chroms=include_chroms)
        return geneinfo
    else:
        try:
            with gzip.open(annotfile, 'r') as mfile:
                for line in mfile:
                    linesplit = line.decode().strip().split('\t')
                    if linesplit[0][0] == '#' or linesplit[2] != feature: continue # skip header

                    chrom = linesplit[0][3:]
                    if include_chrom > 0:
                        include_chroms = ['{:d}'.format(include_chrom)]
                    if chrom not in include_chroms: continue

                    # Any particular biotype selected?
                    infolist = linesplit[8].split(';')

                    if mode == "v19":
                        if len(mybiotype) > 0:
                            rowtype = infolist[2].strip().split(' ')[1].replace('"','')
                            if rowtype not in mybiotype: continue
                        gene_name = infolist[4].strip().split(' ')[1].replace('"','')

                    if mode == "v26":
                        if len(mybiotype) > 0:
                            rowtype = infolist[1].strip().split(' ')[1].replace('"','')
                            if rowtype not in mybiotype: continue
                        gene_name = infolist[2].strip().split(' ')[1].replace('"','')

                    # TSS: gene start (0-based coordinates for BED)
                    if linesplit[6] == '+':
                        start = np.int64(linesplit[3]) - 1
                        end   = np.int64(linesplit[4])
                    elif linesplit[6] == '-':
                        start = np.int64(linesplit[3])  # last base of gene
                        end   = np.int64(linesplit[4]) - 1
                    else:
                        raise ValueError('Strand not specified.')

                    # For simulation
                    if linesplit[1] == 'SIMULATION':
                        start = np.int64(linesplit[3])
                        end   = np.int64(linesplit[4])

                    gene_id = infolist[0].strip().split(' ')[1].replace('"','')
                    if trim:
                        gene_id = gene_id.split(".")[0]
                    this_gene = GeneInfo(name       = gene_name,
                                         ensembl_id = gene_id,
                                         chrom      = int(chrom),
                                         start      = start,
                                         end        = end)

                    geneinfo.append(this_gene)
        except IOError as err:
            raise IOError('{:s}: {:s}'.format(annotfile, err.strerror))

        return geneinfo

def affy_exon_chip(filepath, include_chrom = 0, include_chroms=['{:d}'.format(x + 1) for x in range(22)]):
    geneinfo = list()
    try:
        with open(filepath, 'r') as mfile:
            next(mfile) # skip header
            for line in mfile:
                linesplit = line.strip().split('\t')
                if linesplit[0][0] == '#' : continue 

                chrom = linesplit[2][3:]
                if include_chrom > 0:
                    include_chroms = ['{:d}'.format(include_chrom)]
                if chrom not in include_chroms: continue

                # TSS: gene start (0-based coordinates for BED)
                if linesplit[3] == '+':
                    start = np.int64(linesplit[4]) - 1
                    end   = np.int64(linesplit[5])
                elif linesplit[3] == '-':
                    start = np.int64(linesplit[4])  # last base of gene
                    end   = np.int64(linesplit[5]) - 1
                else:
                    raise ValueError('Strand not specified.')

                gene_name = linesplit[7].split("//")[0].rstrip()
                transcript_cluster_id = linesplit[0]
                this_gene = GeneInfo(name       = gene_name,
                                     ensembl_id = transcript_cluster_id,
                                     chrom      = int(chrom),
                                     start      = start,
                                     end        = end)

                geneinfo.append(this_gene)
    except IOError as err:
        raise IOError('{:s}: {:s}'.format(annotfile, err.strerror))

    return geneinfo