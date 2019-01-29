# Splits a vcf file in chromosomes
# Filters out variants
# Converts snp id to rsids

import os, re
import gzip
import collections
import argparse

def parse_args():

    parser = argparse.ArgumentParser(description='Split vcf file in chromosomes, optional: filter variants and select samples')

    parser.add_argument('--input',
                        type=str,
                        dest='infile',
                        metavar='FILE',
                        help='input VCF file of all chromosomes together')

    parser.add_argument('--rsid-table',
                        type=str,
                        dest='annotfile',
                        metavar='FILE',
                        help='Lookup table file for RSIDs from GTEx')

    parser.add_argument('--outprefix',
                        type=str,
                        dest='outprefix',
                        metavar='FILE',
                        help='output prefix for each chromosome, including output directory (appends .chr__.vcf.gz)')

    parser.add_argument('--filter-snps',
                        dest='filtersnps',
                        action='store_true',
                        help='Filter variants, default False')

    parser.add_argument('--incl-samples',
                        type=str,
                        dest='samplefile',
                        metavar='STR',
                        help='optional: list of samples to keep')

    opts = parser.parse_args()
    return opts



if __name__ == '__main__':
    
    opts = parse_args()

    annotfilepath = opts.annotfile
    infile = opts.infile
    famfile = opts.samplefile
    outfile = opts.outprefix + "chr{:d}.vcf.gz"

    SNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

    rsidlist = collections.defaultdict(lambda:None)
    with open(annotfilepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            linesplit = line.strip().split()
            varid = linesplit[2]
            rsidlist[varid] = linesplit[5]
    print ("Annotation reading complete")

    samples = list()
    if famfile != None:
        with open(famfile) as instream:
            samples = [ line.strip() for line in instream ]

    filter_samples = True if len(samples) else False

    header_lines = list()
    snp_lines = list()
    prev_chrom = None

    counter = 0
    with gzip.open(infile, 'r') as vcf:
        for line in vcf:
            linestrip = line.decode().strip()
            if linestrip[:2] == '##': 
                header_lines.append(line)
                continue
            if linestrip[:6] == '#CHROM':
                linesplit = linestrip.split("\t")
                donor_ids = linesplit[9:]

                newdonors = list()
                #Fix/shorten donor ids
                for d in donor_ids:
                    newdonors.append("-".join(d.strip().split("-")[0:2]))

                # Select samples supplied by user
                if filter_samples:
                    common = [s for s in newdonors if s in samples]
                    if len(common) != len(samples):
                        raise ValueError("Some samples were not found! check the sample file")

                    ix = [newdonors.index(s) for s in samples]
                    sorted_newdonors = [newdonors[i] for i in ix]

                    newdonors_line = "\t".join(sorted_newdonors)
                else:
                    newdonors_line = "\t".join(newdonors)

                new_header_line = "\t".join(linesplit[:9]) + "\t" + newdonors_line + "\n"
                header_lines.append(new_header_line.encode('utf-8'))
            else:
                linesplit = linestrip.split("\t")
                chrom = int(linesplit[0])
                pos   = int(linesplit[1])
                varid = linesplit[2]
                
                # Convert ids to dbSNP142
                rsid  = rsidlist[varid]
                if rsid == None:
                    print("Skipped variant")
                    continue
                ref   = linesplit[3]
                alt   = linesplit[4]
                QC    = linesplit[6]
                QC_meta = linesplit[7].split(";")
                maf = float(QC_meta[0].split("=")[1])

                if prev_chrom == None:
                    prev_chrom = chrom
                    indels_filter=0
                    complement_filter=0
                    pass_filter=0
                    maf_filter=0

                if opts.filtersnps:
                    # Skip indels
                    if len(ref) > 1 or len(alt) > 1:
                        indels_filter +=1
                        continue
                    # Skip ambiguous strands
                    if SNP_COMPLEMENT[ref] == alt:
                        complement_filter +=1
                        continue
                    if QC != "PASS":
                        pass_filter +=1
                        continue
                    if maf < 0.1:
                        maf_filter +=1
                        continue

                if filter_samples:
                    GTs = linesplit[9:]
                    sorted_GTs = [GTs[i] for i in ix]
                    newline = "\t".join([str(chrom), str(pos), rsid] + linesplit[3:9] + sorted_GTs) + "\n"
                else:
                    newline = "\t".join([str(chrom), str(pos), rsid] + linesplit[3:]) + "\n" 

                if chrom == prev_chrom:
                    snp_lines.append(newline.encode('utf-8'))
                else:
                    print("Processed Chromosome {:d}, writing...".format(int(prev_chrom)), end=" ")
                    with gzip.open(outfile.format(int(prev_chrom)), 'wb') as outvcf:
                        for l in header_lines:
                            outvcf.write(l)
                        for s in snp_lines:
                            outvcf.write(s)
                    print("Done")
                    print("{:d} indels deleted".format(indels_filter))
                    print("{:d} complement snps deleted".format(complement_filter))
                    print("{:d} SNPs did not pass QC".format(pass_filter))
                    print("{:d} SNPs with MAF < 0.1".format(maf_filter))
                    snp_lines = list()
                    snp_lines.append(newline.encode('utf-8'))
                    prev_chrom = chrom
                    indels_filter=0
                    complement_filter=0
                    pass_filter=0
                    maf_filter=0
        print("Processed Chromosome {:d}, writing...".format(int(chrom)), end=" ")
        with gzip.open(outfile.format(int(chrom)), 'wb') as outvcf:
            for l in header_lines:
                outvcf.write(l)
            for s in snp_lines:
                outvcf.write(s)
        print("Done")
        print("{:d} indels deleted".format(indels_filter))
        print("{:d} complement snps deleted".format(complement_filter))
        print("{:d} SNPs did not pass QC".format(pass_filter))
        print("{:d} SNPs with MAF < 0.1".format(maf_filter))
        snp_lines = list()
        prev_chrom = chrom