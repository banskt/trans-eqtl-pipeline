import gzip
import argparse
import collections

def parse_args():

    parser = argparse.ArgumentParser(description='Convert vcf.gz files to dosage')

    parser.add_argument('--vcf',
                         help='input file path for vcf.gz',
                         type=str,
                         dest='vcf_file')

    parser.add_argument('--rsid-table',
                         dest='annot_file',
                         type=str,
                         help='annotation file of SNPs')

    parser.add_argument('--out',
                        dest='out_file',
                        help='output file path for dosage')    

    opts = parser.parse_args()
    return opts


if __name__ == '__main__':
    args = parse_args()

    if args.annot_file != None:
        rsidlist = collections.defaultdict(lambda:0)
        with open(args.annot_file, 'r') as mfile:
            next(mfile)
            for line in mfile:
                linesplit = line.strip().split()
                varid = linesplit[2]
                rsidlist[varid] = linesplit[5]
        print ("Annotation reading complete")

    dosagefile = gzip.open(args.out_file, 'w')

    with gzip.open(args.vcf_file, 'r') as vcf:
        for line in vcf:
            if line.decode().strip()[0] == '#': continue
            linesplit = line.decode().strip().split()
            chrom = linesplit[0]
            varid = linesplit[2]
            rsid = rsidlist[varid] if args.annot_file != None else varid
            pos = linesplit[1]
            ref = linesplit[3]
            alt = linesplit[4]
            snpdosage = [float(x.split(':')[2]) for x in linesplit[9:]]
            freq = sum(snpdosage) / 2 / len(linesplit[9:])
            maf = freq
            # Only used for GEUVADIs study (on predixcan)
            # if freq > 0.5:
            #     maf = 1 - freq
            #     ref = linesplit[4]
            #     alt = linesplit[3]
            #     snpdosage = [2 - x for x in snpdosage]
            snpdosage = ['{:g}'.format(x) for x in snpdosage]
            dosagefile.write('{:s} {:s} {:s} {:s} {:s} {:g} {:s}\n'.format(chrom, rsid, pos, ref, alt, maf, ' '.join(snpdosage)).encode())

    dosagefile.close()
