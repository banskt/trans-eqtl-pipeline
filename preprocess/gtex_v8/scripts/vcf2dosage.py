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

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    
    opts = parse_args()
    vcffile = opts.invcf

    SNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

    dosage = list()
    snpinfo = list()
    linenum = 0
    mode = "DS"
    with gzip.open(opts.outvcf, 'wb') as outvcf:
        with gzip.open(vcffile, 'r') as vcf:
            for line in vcf:
                linestrip = line.decode().strip()
                if linestrip[:2] == '##': 
                    continue
                    # outvcf.write(line)
                if linestrip[:6] == '#CHROM':
                    continue
                    # linesplit = linestrip.split("\t")
                    # donor_ids = linesplit[9:]
                    # outvcf.write(line)
                else:
                    linesplit = linestrip.split("\t")
                    if linesplit[0].startswith("chr"):
                        chrom = int(linesplit[0][3:])
                    else:
                        chrom = int(linesplit[0])
                    pos   = int(linesplit[1])
                    varid = linesplit[2]
                    ref   = linesplit[3]
                    alt   = linesplit[4]

                    if mode == "DS":
                        if "DS" not in linesplit[8].split(':'):
                            mode = "GT"
                        else:
                            dsindx = linesplit[8].split(':').index("DS")
                            ds = [x.split(':')[dsindx] if len(x) > 1 else "." for x in linesplit[9:]]
                            gtindx = linesplit[8].split(':').index("GT")
                            for i, x in enumerate(ds):
                                if x == ".":
                                    gt = linesplit[9+i].split(':')[gtindx]
                                    if len(gt) == 3 and gt[0] != "." and gt[2] != ".":
                                        ds[i] = float(int(gt[0]) + int(gt[2]))

                    if mode == "GT":
                        if "GT" not in linesplit[8].split(':'):
                            print("ERROR: no GT field in VCF file")
                            raise
                        gtindx = linesplit[8].split(':').index("GT")
                        gt = [x.split(':')[gtindx] if len(x) > 1 else "." for x in linesplit[9:]]
                        ds = [ float(int(x[0]) + int(x[2])) if len(x) == 3 and x[0] != "." and x[2] != "." else "." for x in gt ]

                    ds_notna = [float(x) for x in ds if x != "."]
                    freq = sum(ds_notna) / 2 / len(ds_notna)
                    maf = freq
                    snpdosage = [float(x) if x != '.' else 2 * freq for x in ds]

                    linenum+=1
                    # if linenum > 5000:
                    #     break
                    outvcf.write('{:d} {:s} {:d} {:s} {:s} {:g} {:s}\n'.format(chrom, varid, pos, ref, alt, maf, ' '.join(snpdosage)).encode())