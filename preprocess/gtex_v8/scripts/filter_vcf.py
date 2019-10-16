import os, re
import gzip
import collections
import argparse


def parse_args():

    parser = argparse.ArgumentParser(description='Filter VCF file.')

    parser.add_argument('--in',
                        dest='invcf',
                        metavar='FILE',
                        help='input vcf')

    parser.add_argument('--out',
                        dest='outvcf',
                        metavar='FILE',
                        help='output vcf')

    parser.add_argument('--maf',
                    type=float,
                    dest='maf_thres',
                    metavar='FLOAT',
                    help='output vcf')

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    
    opts = parse_args()
    vcffile = opts.invcf
    maf_cutoff=opts.maf_thres

    SNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

    dosage = list()
    snpinfo = list()
    linenum = 0
    indels_filter=0
    complement_filter=0
    maf_filter=0
    with gzip.open(opts.outvcf, 'wb') as outvcf:
        with gzip.open(vcffile, 'r') as vcf:
            for line in vcf:
                linestrip = line.decode().strip()
                if linestrip[:2] == '##': 
                    outvcf.write(line)
                    continue
                if linestrip[:6] == '#CHROM':
                    linesplit = linestrip.split("\t")
                    donor_ids = linesplit[9:]
                    outvcf.write(line)
                else:
                    linesplit = linestrip.split("\t")
                    chrom = int(linesplit[0][3:])
                    pos   = int(linesplit[1])
                    varid = linesplit[2]
                    ref   = linesplit[3]
                    alt   = linesplit[4]

                    gtindx = linesplit[8].split(':').index("GT")
                    gt = [x.split(':')[gtindx] for x in linesplit[9:]]
                    ds = [ float(int(x[0]) + int(x[2])) if len(x) == 3 and x[0] != "." and x[2] != "." else "." for x in gt ]

                    ds_notna = [float(x) for x in ds if x != "."]
                    freq = sum(ds_notna) / 2 / len(ds_notna)
                    maf = freq
                    snpdosage = [float(x) if x != '.' else 2 * freq for x in ds]

                    linenum+=1
                    # if linenum > 5000:
                    #     break
                    # Skip indels
                    if len(ref) > 1 or len(alt) > 1:
                        indels_filter +=1
                        continue
                    # Skip ambiguous strands
                    if SNP_COMPLEMENT[ref] == alt:
                        complement_filter +=1
                        continue
                    if maf < maf_cutoff or maf > (1 - maf_cutoff):
                        maf_filter +=1
                        continue

                    outvcf.write(line)
                    
        print("{:d} indels deleted".format(indels_filter))
        print("{:d} complement snps deleted".format(complement_filter))
        print("{:d} SNPs with MAF < {:d}".format(maf_filter, maf_cutoff))