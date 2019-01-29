# select samples from a single vcf file
# Author: Franco L. Simonetti

import os
import re
import gzip
import collections
import argparse

def parse_args():

    parser = argparse.ArgumentParser(description='Select and sort samples for a given VCF chrom file')

    parser.add_argument('--input',
                        type=str,
                        dest='infile',
                        metavar='FILE',
                        help='input VCF file of a single chromosome')

    parser.add_argument('--out',
                        type=str,
                        dest='outfile',
                        metavar='FILE',
                        help='output VCF file')

    parser.add_argument('--incl-samples',
                        type=str,
                        dest='samplefile',
                        metavar='STR',
                        help='list of samples to keep and sorting order')

    parser.add_argument('--annot',
                        type=str,
                        dest='annotfile',
                        metavar='STR',
                        help='SNP annotation file with two tab-separated columns -- original annotation and new annotation')

    opts = parser.parse_args()
    return opts


if __name__ == '__main__':
    
    opts = parse_args()

    infile = opts.infile
    famfile = opts.samplefile
    outfile = opts.outfile

    with open(famfile, 'r') as instream:
        gx_samples = [ line.strip() for line in instream ]

    change_annot = False
    if opts.annotfile is not None:
        change_annot = True
        rsid_annot = dict()
        with open(opts.annotfile, 'r') as instream:
            for line in instream:
                linesplit = line.strip().split()
                rsid_annot[linesplit[0]] = linesplit[1]

    header_lines = list()
    snp_lines = list()
    n_no_lookup = 0

    with gzip.open(infile, 'r') as vcfstream:
        for line in vcfstream:
            linestrip = line.decode().strip()
            if linestrip[:2] == '##': 
                header_lines.append(line)
                continue
            if linestrip[:6] == '#CHROM':
                linesplit = linestrip.split("\t")
                vcf_samples = linesplit[9:]

                # Select samples supplied by user
                common_samples = [s for s in vcf_samples if s in gx_samples]
                if len(common_samples) != len(gx_samples):
                    raise ValueError("Some samples were not found! check the sample file")
                ix_sorted = [vcf_samples.index(s) for s in gx_samples]
                vcf_samples_keep = [vcf_samples[i] for i in ix_sorted]

                thisline = "\t".join(linesplit[:9]) + "\t" + "\t".join(vcf_samples_keep) + "\n"
                #new_header_line = "\t".join(linesplit[:9]) + "\t" + newdonors_line + "\n"
                header_lines.append(thisline.encode('utf-8'))
            else:
                linesplit = linestrip.split("\t")
                if change_annot:
                    snpid = linesplit[2]
                    try:
                        linesplit[2] = rsid_annot[snpid]
                    except:
                        n_no_lookup += 1
                gt = linesplit[9:]
                gt_sorted = [gt[i] for i in ix_sorted]
                newline = "\t".join(linesplit[:9]) + "\t" + "\t".join(gt_sorted) + "\n"
                snp_lines.append(newline.encode('utf-8'))

    outstring = f"Found {len(snp_lines):d} SNPs for {len(common_samples):d} samples"
    if change_annot:
        print( outstring + f", no matching rsid for {n_no_lookup:d} SNPs" )
    else:
        print( outstring )

    with gzip.open(outfile, 'wb') as outvcf:
        for l in header_lines:
            outvcf.write(l)
        for s in snp_lines:
            outvcf.write(s)
