import os
import argparse 

def parse_args():

    parser = argparse.ArgumentParser(description='Add SNP position to MatrixEQTL result files')

    parser.add_argument('--input',
                        dest='infile',
                        metavar='FILE',
                        help='directory of MatrixEQTL chromosome result files (trans_eqtls.txt and cis_eqtls.txt')

    parser.add_argument('--posfile',
                        dest='posfile',
                        metavar="FILE",
                        help='SNP pos file for all chromosomes')


    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    opts = parse_args()

    # load SNP rsid position dictionary
    snpdict = dict()
    with open(opts.posfile) as instream:
        for line in instream:
            arr = line.rstrip().split()
            snpdict[arr[1]] = arr[2]

    with open(opts.infile+".pos", 'w') as outstream:
        with open(opts.infile) as instream:
            next(instream)
            for line in instream:
                newline = line.rstrip()
                rsid    = newline.split()[0]
                if snpdict.get(rsid, False):
                    outstream.write(newline+"\t"+snpdict[rsid]+"\n")
                else:
                    print("ERROR: RSID not found")
                    raise

    os.remove(opts.infile)
        