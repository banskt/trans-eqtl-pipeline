# select samples from a single vcf file
# Author: Franco L. Simonetti

import gzip
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

    opts = parser.parse_args()
    return opts


if __name__ == '__main__':
    
    opts = parse_args()

    infile = opts.infile
    famfile = opts.samplefile
    outfile = opts.outfile

    with open(famfile, 'r') as instream:
        gx_samples = [ line.strip() for line in instream ]


    nsnps = 0
    with gzip.open(infile, 'r') as vcfstream, gzip.open(outfile, 'wb') as outstream:
        for line in vcfstream:
            linestrip = line.decode().strip()
            if linestrip[:2] == '##': 
                outstream.write(line)
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
                outstream.write(thisline.encode('utf-8'))
            else:
                linesplit = linestrip.split("\t")
                gt = linesplit[9:]
                gt_sorted = [gt[i] for i in ix_sorted]
                newline = "\t".join(linesplit[:9]) + "\t" + "\t".join(gt_sorted) + "\n"
                outstream.write(newline.encode('utf-8'))
                nsnps += 1

    outstring = f"Found {nsnps:d} SNPs for {len(common_samples):d} samples"
    print (outstring)
