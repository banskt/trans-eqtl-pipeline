import argparse
import gzip
import numpy as np

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

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':

    opts = parse_args()
    infile = opts.infile
    outfile = opts.outfile
 
    GT = list(['0/0', '0/1', '1/1'])

    with gzip.open(opts.infile, 'r') as instream, gzip.open(opts.outfile, 'w') as outstream:
        for line in instream:
            linestrip = line.decode().strip()
            if linestrip[0] == '#':
                outstream.write(line)
            else:
                linesplit = linestrip.split("\t")
                genotypes = linesplit[9:]
                gt = [g.split(':') for g in genotypes]
                for i, g in enumerate(gt):
                    if g[0] == './.':
                        #gt[i][0] = GT[np.argmax(list(map(float, g[1].split(","))))]
                        gt[i][0] = GT[np.argmax(np.array(g[1].split(",")).astype(float))]
                linesplit[9:] = [":".join(g) for g in gt]
                newline = "\t".join(linesplit) + "\n"
                outstream.write(newline.encode('utf-8'))
