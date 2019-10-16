import gzip
import numpy as np
import re
import argparse

def parse_args():

    parser = argparse.ArgumentParser(description='Collaspe phASER RNA-seq phased counts')

    parser.add_argument('--in',
                        type=str,
                        dest='infile',
                        metavar='FILE',
                        help='Input allele-specific phASER matrix file')

    parser.add_argument('--out',
                        type=str,
                        dest='outfile',
                        metavar='FILE',
                        help='Output collapsed counts file')


    opts = parser.parse_args()
    return opts


opts = parse_args()

inmatx = opts.infile #"/cbscratch/franco/datasets/gtex_v8/expression/phASER_GTEx_v8_matrix.txt.gz"
outmatx = opts.outfile #"/cbscratch/franco/datasets/gtex_v8/expression/phASER_GTEx_v8_matrix_collapsed_counts_new.txt.gz"

l=0
with gzip.open(outmatx, 'wb') as outstream:
    with gzip.open(inmatx) as instream:
        header = instream.readline().decode()
        harr   = header.split()
        sample_ids = harr[4:]
        idline = "\t".join(sample_ids)
        hline  = harr[1]+"\t"+idline+"\n"
        outstream.write(hline.encode())
        for line in instream:
            arr = line.decode().rstrip().split()
            if re.search(r'[XYM]', arr[0][3:]):
                continue
            counts = [np.sum(np.array(x.split("|"), dtype=int)) for x in arr[4:] ]
            countline = arr[1]+"\t"+"\t".join([str(i) for i in counts])+"\n"
            outstream.write(countline.encode())
            l += 1
            if (l % 1000) == 0:
                print("Processed {:d} genes".format(l))