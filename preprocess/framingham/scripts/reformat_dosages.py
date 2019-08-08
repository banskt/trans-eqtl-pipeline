import gzip, os, sys
import collections
import pandas as pd
import argparse
import re
import time
import numpy as np

BASES = ["A", "C", "G", "T"]
SNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

SNPINFO_FIELDS = ['chrom', 'varid', 'bp_pos', 'ref_allele', 'alt_allele', 'maf']
class SnpInfo(collections.namedtuple('_SnpInfo', SNPINFO_FIELDS)):
    __slots__ = ()

def parse_args():
    parser = argparse.ArgumentParser(description='Reformats Framingham genotype files')

    parser.add_argument('--gt',
                         help='.gz genotype file of dosages (one line per sample, one column per SNP)',
                         type=str,
                         dest='gt_file')

    parser.add_argument('--info',
                         dest='info_file',
                         type=str,
                         help='.gz file with SNP info')

    parser.add_argument('--sample',
                         dest='sample_file',
                         type=str,
                         help='file with sample info manifest') #"../phg000679.v2_release_manifest.txt"

    parser.add_argument('--out',
                        dest='out_file',
                        default=None,
                        help='output file for genotype')    

    parser.add_argument('--consent',
                        dest='c',
                        type=str,
                        default='c2',
                        help='consent group (either c1 or c2)')

    parser.add_argument('--keepfile',
                        dest='keepfile',
                        type=str,
                        default=None,
                        help='File with snps retained in Joehannes et al') 

    parser.add_argument('--chrm',
                        dest='chrom',
                        type=int,
                        default=None,
                        help='Chromosome number') 

    opts = parser.parse_args()
    return opts

def read_keepfile(keepfile, chrm):
    posdict = collections.defaultdict(lambda: False)
    dups = 0
    with open(keepfile) as instream:
        next(instream)
        for line in instream:
            linestr = line.rstrip()
            arr = linestr.split(":")
            this_chrm = int(arr[0]) if arr[0] != 'X' else 999
            if this_chrm < chrm:
                continue
            if this_chrm == chrm:
                if posdict[linestr]:
                    dups += 1
                posdict[linestr] = True
            if this_chrm > chrm:
                break
    print("Found",dups,"duplicate SNP pos")
    return posdict


def write_dosages(snp_info, dosages, outfile):
    with gzip.open(outfile, 'wb') as outstream:
        for i,snp in enumerate(snp_info):
            # print(i, snp)
            dosage_row = " ".join([str(x) for x in dosages[i,:]])
            dosage_maf = (1 - sum(dosages[i,:]) / 2 / len(dosages[i,:]))
            outstream.write("{:d} {:s} {:d} {:s} {:s} {:g} {:s}\n".format(snp.chrom, snp.varid, snp.bp_pos, snp.ref_allele, snp.alt_allele, dosage_maf, dosage_row).encode('utf-8'))

if __name__ == '__main__':
    args = parse_args()

    sample_df = pd.read_csv(args.sample_file, sep="\t", comment = "#")
    samples = sample_df[sample_df["File_Name_Mtrx"].str.contains(args.c)]

    keepflag = False
    if args.keepfile != None and os.path.exists(args.keepfile):
        posdict = read_keepfile(args.keepfile, args.chrom)
        keepflag = True
    else:
        posdict = collections.defaultdict(lambda: False)
        print("No KEEP file found.")

    i = -1
    nmaf = 0
    nbase = 0
    nambi = 0
    total = 0
    nrsq = 0
    good_snps = []
    good_ids = np.array([], dtype=int)
    keep_snps = []
    keep_ids = np.array([], dtype=int)
    with gzip.open(args.info_file, 'r') as instream:
        next(instream)
        for line in instream:
            i += 1
            arr = line.decode("utf-8").split()
            this_snp = SnpInfo(chrom      = int(arr[0].split(":")[0]),
                               bp_pos     = int(arr[0].split(":")[1]),
                               varid      = "{:s}_{:s}_{:s}".format(arr[0],arr[1], arr[2]),
                               ref_allele = arr[1],
                               alt_allele = arr[2],
                               maf        = float(arr[4]))
            # removes badly imputed snps
            if posdict[arr[0]]:
                keep_snps.append(this_snp)
                keep_ids = np.append(keep_ids, i)
            if float(arr[6]) < 0.3:
                nrsq += 1
                continue
            if this_snp.maf < 0.01 or this_snp.maf > 0.99:
                nmaf += 1
                continue
            if this_snp.ref_allele not in BASES or this_snp.alt_allele not in BASES:
                nbase += 1
                continue
            if SNP_COMPLEMENT[this_snp.ref_allele] == this_snp.alt_allele:
                nambi += 1
                continue
            good_snps.append(this_snp)
            good_ids = np.append(good_ids, i)
            total += 1
    print("low imputation Rsq:", nrsq)
    print("low maf:", nmaf)
    print("not bases:", nbase)
    print("ambiguous:", nambi)
    print("SNPs left:", total)
    print("Keep SNPs:", len(keep_snps))

    s = time.time()

    gt_sample_ids = list() #list(samples["Subject_ID"])
    dosages_arr = np.zeros((samples.shape[0], total))
    if keepflag:
        dosages_keep = np.zeros((samples.shape[0], len(keep_snps)))
    lines = 0
    with gzip.open(args.gt_file) as instream:
        for line in instream:
            arr = line.decode('utf-8').split("\t")
            gt_sample_ids.append(str(arr[0]))
            dosages_arr[lines, :] = np.array(arr[2:], dtype=float)[good_ids]
            if keepflag:
                dosages_keep[lines, :] = np.array(arr[2:], dtype=float)[keep_ids]
            lines += 1
        print("Dosage reading took", time.time() - s, "seconds")

    samples_list = list(sample_df["Subject_ID"])
    ix = np.array([gt_sample_ids.index(str(i)) for i in samples_list if str(i) in gt_sample_ids])
    sorted_gt_sample_ids = [gt_sample_ids[i] for i in ix]
    for i,j in zip(sorted_gt_sample_ids, list(samples["Subject_ID"])):
        if int(i) != j:
            print("Sample error: samples are not same", i, j)
            raise

    s = time.time()
    write_dosages(good_snps, dosages_arr[ix, :].T, args.out_file+".txt.gz")
    print("File written to", args.out_file)
    print("Writting took", time.time()-s)

    if keepflag:
        write_dosages(keep_snps, dosages_keep[ix, :].T, args.out_file+".joehannes.txt.gz")
        print("File written to", args.out_file)
        print("Writting took", time.time()-s)

    with open("{:s}.sample".format(args.out_file), 'w') as outstream:
        outstream.write("ID_1 ID_2 missing father mother pheno\n")
        outstream.write("0 0 0 D D B\n")
        for s in sorted_gt_sample_ids:
            outstream.write("{:s} {:s} 0 0 0 -9\n".format(s, s))

    ### Not same order as GENOTYPE!
    # # Write sample file, only required once for each consent group
    # with open(args.out_file+".fhs.sample", 'w') as outstream:
    #     outstream.write("ID_1 ID_2 missing father mother pheno\n")
    #     outstream.write("0 0 0 D D B\n")
    #     for s in sample_ids:
    #         outstream.write("{:d} {:d} 0 0 0 -9\n".format(s, s))
