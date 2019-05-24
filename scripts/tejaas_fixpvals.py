import numpy as np
from scipy import stats
import argparse
import shutil

def parse_args():

    parser = argparse.ArgumentParser(description = "Fix the p-values from TEJAAS")

    parser.add_argument('--infile',
                        type = str,
                        dest = 'infile',
                        metavar = 'FILE',
                        help = 'old TEJAAS result file')

    parser.add_argument('--outfile',
                        type = str,
                        dest = 'outfile',
                        metavar = 'FILE',
                        help = 'corrected TEJAAS result file')

    opts = parser.parse_args()
    return opts

def read_tejaas(filename):
    rsidlist = list()
    bpposlist = list()
    qscrlist = list()
    qmeanlist = list()
    qvarlist  = list()
    qscalelist = list()
    with open(filename, 'r') as instream:
        header = next(instream)
        for line in instream:
            linesplit = line.strip().split()
            rsid = linesplit[0]
            bppos = linesplit[1]
            qscr = float(linesplit[2])
            qmean = float(linesplit[3])
            qvar = float(linesplit[4])
            qscale = (qscr - qmean) / qvar
            rsidlist.append(rsid)
            bpposlist.append(bppos)
            qscrlist.append(qscr)
            qmeanlist.append(qmean)
            qvarlist.append(qvar)
            qscalelist.append(qscale)
    return header, rsidlist, bpposlist, qscrlist, qmeanlist, qvarlist, qscalelist

if __name__ == '__main__':

    opts = parse_args()

    infile = opts.infile
    outfile = opts.outfile

    if infile == outfile:
        backupfile = "{:s}.backup".format(opts.infile)
        shutil.copy (infile, backupfile)

    header, rsids, bppos, qstats, qmeans, qvars, qscales = read_tejaas(infile)
    pvals = 2.0 * (1 - stats.norm.cdf(np.abs(qscales)))
    with open(outfile, 'w') as fout:
        fout.write(header)
        for i, rsid in enumerate(rsids):
            line = "{:s}\t{:s}\t{:g}\t{:g}\t{:g}\t{:g}\n".format(rsid, bppos[i], qstats[i], qmeans[i], qvars[i], pvals[i])
            fout.write(line)
