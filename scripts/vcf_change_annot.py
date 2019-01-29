import os
import gzip
import argparse


def parse_args():

    parser = argparse.ArgumentParser(description='Select and sort samples for a given VCF chrom file')

    parser.add_argument('--input',
                        type=str,
                        dest='infile',
                        metavar='FILE',
                        help='input gzipped VCF file of a single chromosome')

    parser.add_argument('--out',
                        type=str,
                        dest='outfile',
                        metavar='FILE',
                        help='output gzipped VCF file')

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
    outfile = opts.outfile
    annotfile = opts.annotfile

    rsid_annot = dict()
    with open(opts.annotfile, 'r') as instream:
        for line in instream:
            linesplit = line.strip().split()
            rsid_annot[linesplit[0]] = linesplit[1]

    n_no_lookup = 0
    nsnps = 0
    with gzip.open(infile, 'r') as vcfstream, gzip.open(outfile, 'w') as outstream:
        for line in vcfstream:
            linestrip = line.decode().strip()
            if linestrip[0] == '#':
                outstream.write(line)
            else:
                linesplit = linestrip.split("\t")
                nsnps += 1
                snpid = linesplit[2]
                try:
                    linesplit[2] = rsid_annot[snpid]
                except:
                    n_no_lookup += 1
                newline = "\t".join(linesplit) + "\n"
                outstream.write(newline.encode('utf-8'))

    outstring = f"{nsnps:d} SNPs after filtering, no matching rsid for {n_no_lookup:d} SNPs"
    print(outstring)
