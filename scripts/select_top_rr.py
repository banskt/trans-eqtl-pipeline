import argparse
import gzip
import pandas as pd
import os 

def parse_args():

    parser = argparse.ArgumentParser(description='Filter top trans-eQTLs from RR output.')

    parser.add_argument('--input',
                        type=str,
                        dest='infile',
                        metavar='FILE',
                        help='rr.txt compiled chunks file')

    parser.add_argument('--n',
                        type=int,
                        dest='nsnps',
                        default=1000,
                        help='number N of best snps to keep')


    opts = parser.parse_args()
    return opts

if __name__ == '__main__':

    opts = parse_args()

    if os.path.exists(opts.infile):
        df = pd.read_csv(opts.infile, header=0, sep="\t")
        outdf = df.sort_values(by=['P']).iloc[:opts.nsnps].sort_values(by=['Pos'])
        with open(opts.infile+".top"+str(opts.nsnps), 'w') as outstream:
            outdf.to_csv(outstream, sep="\t", header=True, index=False)
    else:
        print(opts.infile, "does not exist")
        raise