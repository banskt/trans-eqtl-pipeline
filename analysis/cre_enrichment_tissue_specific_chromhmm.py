import os
import numpy as np
import argparse

from utils import utils
from utils import read_tejaas_results

CHROMATIN_REGIONS = ['Promoter', 'Enhancer', 'Transcribed', 'ZNF', 'Heterochromatin', 'Inactive', 'Quiescent', 'LowConfidence']
CHROMATIN_REGIONS = ['Active', 'Inactive', 'Silent']


def get_label_dict():
    ## returns the index of the global CHROMATIN_REGIONS list
    labels = dict() 
    labels['1'] = 0
    labels['2'] = 0
    labels['3'] = 2
    labels['4'] = 2
    labels['5'] = 2
    labels['6'] = 1
    labels['7'] = 1
    labels['8'] = 3
    labels['9'] = 4
    labels['10'] = 5
    labels['11'] = 5
    labels['12'] = 5
    labels['13'] = 5
    labels['14'] = 5
    labels['15'] = 6
    return labels


def parse_args():

    parser = argparse.ArgumentParser(description='Calculate tissue-matched enrichment of cis-regulatory elements')

    parser.add_argument('--method',
                        type=str,
                        dest='method',
                        default='chromhmm',
                        help='name of segmentation method (segway or chromhmm)')

    parser.add_argument('--dhsdir',
                        type=str,
                        dest='dhsdir',
                        metavar='FILE',
                        help='full path of the directory containing DHS regions')

    parser.add_argument('--tissue',
                        type=str,
                        dest='tissuefile',
                        metavar='FILE',
                        help='full path of the tissue file')

    parser.add_argument('--out',
                        type=str,
                        dest='outfile',
                        metavar='FILE',
                        help='name of the output file')

    parser.add_argument('--resdir',
                        type=str,
                        dest='resdir',
                        metavar='DIR',
                        help='name of the result directory')

#    parser.add_argument('--tissue',
#                        nargs='*',
#                        type=str,
#                        dest='tissuelist',
#                        metavar='TISSUE',
#                        help='short name of all tissues to be analysed')


    parser.add_argument('--match',
                        type=str,
                        dest='matchfile',
                        metavar='FILE',
                        help='name of the tissue matching file')
 

    parser.add_argument('--eidfile',
                        type=str, 
                        dest='eidfile', 
                        metavar='FILE',
                        help='name of the file containing epigenome IDs and their mapping to DHS headers')

    opts = parser.parse_args()
    return opts


def get_epigenome_dict(filename, method):
    epigenomes = dict()
    with open(filename, 'r') as instream:
        for line in instream:
            lsplit = line.split(";")
            if method == 'segway':
                indx = 2
            elif method == 'chromhmm':
                indx = 1
            eid = lsplit[0].rstrip()
            epigenomes[eid] = lsplit[indx].rstrip()
    return epigenomes


def get_random_resdicts(dirname):
    rand_dicts = list()
    for i in range(10):
        res_dict = dict()
        iterdir = os.path.join(dirname, f'random_50000_{i+1 :02d}')
        for chrm in range(1, 23):
            chrm_bppos_list = list()
            filename = os.path.join(iterdir, f'chr{chrm}.txt')
            with open(filename, 'r') as infile:
                for line in infile:
                    line = line.strip()
                    spos = int(line.split()[0].strip())
                    chrm_bppos_list.append(spos)
            chrm_bppos_list.sort()
            res_dict[chrm] = chrm_bppos_list
        rand_dicts.append(res_dict)
    return rand_dicts


def find_dhs_overlap(res_dict, dhs_file, outfile):
    dhs = open(dhs_file)
    line = dhs.readline()
    prev_chrm = 0
    nannot = 0
    while line:
        arr = line.rstrip().split()
        chrm = int(arr[0][3:])
        start = int(arr[1])
        end = int(arr[2])
        if chrm != prev_chrm:
            remaining = res_dict[chrm]
            checked = 0
        if len(remaining) == 0:
            ## No more SNPs in this chromosome, just continue reading the DHS file
            line = dhs.readline()
        else:
            for pos in remaining:
                if pos < start:
                    checked += 1
                    remaining = res_dict[chrm][checked:]
                    continue # go to next SNP
                elif pos > end:
                    line = dhs.readline()
                    break # go to next DHS line, keep checking the remaining results
                else:
                    # this is an annotated SNP
                    checked += 1
                    remaining = res_dict[chrm][checked:]
                    nannot += 1
                    outfile.write(line)
                    continue # go to next SNP
        prev_chrm = chrm
    dhs.close()
    return nannot


def read_dhsoverlap_outfile(overlapfilelist, labeldict):
    tannot = np.zeros(len(CHROMATIN_REGIONS))
    for overlapfile in overlapfilelist:
        with open(overlapfile, 'r') as infile:
            for line in infile:
                lsplit = line.split("\t")
                chrom_mark = lsplit[3].strip()
                tannot[labeldict[chrom_mark]] += 1
    return tannot


def get_chromatin_overlap(res_dict, rand_dicts, dhsdir, eidlist, labeldict, tmpoutdir):

    filelist = list()
    for eid in eidlist:
        tmp_overlap_file = os.path.join(tmpoutdir, f'overlapped_regions_{eid}.txt')
        #dhsfile = os.path.join(dhsdir, f'{eid}_15_coreMarks_hg38lift_stateno_clean.bed')
        #ofile = open(tmp_overlap_file, 'w')
        #nannot = find_dhs_overlap(res_dict, dhsfile, ofile)
        #ofile.close()
        filelist.append(tmp_overlap_file)

    tannot = read_dhsoverlap_outfile(filelist, labeldict)
    #tannot = np.zeros(len(CHROMATIN_REGIONS))

    tannot_rand = [np.zeros(len(CHROMATIN_REGIONS)) for x in rand_dicts]
    for i,rdict in enumerate(rand_dicts):
        filelist = list()
        for eid in eidlist:
            tmp_overlap_file = os.path.join(tmpoutdir, f'rand_overlapped_regions_{i+1 :02d}_{eid}.txt')
            #dhsfile = os.path.join(dhsdir, f'{eid}_15_coreMarks_hg38lift_stateno_clean.bed')
            #ofile = open(tmp_overlap_file, 'w')
            #nannot = find_dhs_overlap(rdict, dhsfile, ofile)
            #ofile.close()
            filelist.append(tmp_overlap_file)

        tannot_rand[i] = read_dhsoverlap_outfile(filelist, labeldict)

    return tannot, tannot_rand
    

random_snp_dir = "/usr/users/sbanerj/gtex_v8/genotype/all_samples/random_sampling"

if __name__ == '__main__':
    opts = parse_args()
#    tissuelist = opts.tissuelist
    tissuefile = opts.tissuefile
    resdir = opts.resdir
    dhsdir = opts.dhsdir
    matchfile = opts.matchfile
    eidfile = opts.eidfile
    method = opts.method
    outfile = opts.outfile

    matching_eid  = utils.read_matching_eid(matchfile)
    epigenomedict = get_epigenome_dict(eidfile, method)
    labeldict = get_label_dict()
    tissuelist, tfulls = utils.read_tissues(tissuefile)
    tissue_names = dict()
    for tshort, tfull in zip(tissuelist, tfulls):
        tissue_names[tshort] = tfull

    rand_dicts = get_random_resdicts(random_snp_dir)

    fout = open(outfile, 'w')
    enrichment_string = '\t'.join([x.upper() for x in CHROMATIN_REGIONS])
    fout.write(f'TISSUE\tN_TRANSEQTLS\t{enrichment_string}\n')

    for tissue in tissuelist:
        resfilename = os.path.join(resdir, tissue, 'trans_eqtls.txt')
        transeqtls = read_tejaas_results.transeqtls(resfilename)
        nteqtl = len(transeqtls)
        eidlist = matching_eid[tissue]

        if nteqtl > 0 and len(eidlist) > 0:
            epigenomelist = [epigenomedict[x] for x in eidlist]

            print (f'{tissue_names[tissue]}: {"; ".join(epigenomelist)}')

            res_dict = dict()

            for chrm in range(1, 23):
                res_dict[chrm] = list()

            for teqtl in transeqtls:
                chrm = teqtl.chrom
                spos = teqtl.bp_pos
                res_dict[chrm].append(spos)

            tmpoutdir = os.path.join(os.path.dirname(os.path.abspath(outfile)), tissue)
            if not os.path.exists(tmpoutdir): os.makedirs(tmpoutdir)

            tannot, tannot_rand = get_chromatin_overlap(res_dict, rand_dicts, dhsdir, eidlist, labeldict, tmpoutdir)

            tfrac = tannot / nteqtl
            nx = [sum([len(x) for key, x in rdict.items()]) for rdict in rand_dicts]
            tfrac_rand = [x / nx[i] for i, x in enumerate(tannot_rand)]
            bgfrac = np.mean(np.array(tfrac_rand), axis = 0)
            enrichment = [x / y if y > 0 else np.nan for x, y in zip(tfrac, bgfrac)]
            enrichment_string = '\t'.join([f'{x :g}' if not np.isnan(x) else 'NA' for x in enrichment])
            fout.write(f'{tissue}\t{nteqtl}\t{enrichment_string}\n')
            #print(tissue, enrichment_string[2], enrichment_string[3])
            print(tissue, tannot)
