import time
import sys, os, collections
import numpy as np
from operator import attrgetter
import argparse

def parse_args():

    parser = argparse.ArgumentParser(description='Filter SNPs in LD.')

    parser.add_argument('--infiles',
                        nargs='+',
                        dest='infiles',
                        metavar='FILE',
                        help='Genome-wide trans-eqtls results (all CHRMs in one file - RR or mEQTL format with position)')

    parser.add_argument('--ldfile',
                        type=str,
                        dest='ldfile',
                        help='LD map base filename (eg. chr[CHRM].ld_geno)')

    parser.add_argument('--outfiles',
                        nargs='+',
                        dest='outfiles',
                        help='List of output files corresponding to the input files')


    opts = parser.parse_args()
    return opts


def read_ldfile(ldfile):
    ldict = collections.defaultdict(lambda: False)
    with open(ldfile) as instream:
        next(instream)
        for line in instream:
            arr = line.rstrip().split()
            chrm = arr[0]
            if chrm.startswith("chr"):
                chrm = int(chrm[3:])
            else:
                chrm = int(chrm)
            pos1 = str(arr[1])
            pos2 = str(arr[2])
            n = int(arr[3])
            r2 = float(arr[4])
            if ldict[pos1]:
                ldict[pos1][pos2] = r2
            else:
                ldict[pos1] = collections.defaultdict(lambda: False)
                ldict[pos1][pos2] = r2
                
            if ldict[pos2]:
                ldict[pos2][pos1] = r2
            else:
                ldict[pos2] = collections.defaultdict(lambda: False)
                ldict[pos2][pos1] = r2
    return ldict

import mpmath
mpmath.mp.dps = 500
def pval(x): return float(1 - 0.5 * (1 + mpmath.erf(x/mpmath.sqrt(2))))
   
SNPRES_FIELDS = ['rsid', 'chrom', 'pos', 'q', 'mu', 'sigma', 'p', 'logp', 'maf']
class SNPRes(collections.namedtuple('_SNPRes', SNPRES_FIELDS)):
    __slots__ = ()
        
def tejaas_chr(filepath):
    res = list()
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr   = line.strip().split("\t")
            rsid  = arr[0]
            chrom = int(arr[1])
            pos   = int(arr[2])
            maf   = float(arr[3])
            q     = float(arr[4])
            mu    = float(arr[5])
            sigma = float(arr[6])
            p     = float(arr[7])
            if sigma == 0:
                continue
            p = p if p != 0 else pval( (q - mu) / sigma)
            logp = np.log10(p) # mpmath.log10
            res.append(SNPRes(rsid=rsid, chrom=chrom, pos=pos, q=q, mu=mu, sigma=sigma, p=p, logp=-logp, maf=maf))
    return res

def tejaas_write(snplist, filepath):
    with open(filepath, 'w') as mfile:
        mfile.write("ID\tChrom\tPos\tMAF\tQ\tMu_Q\tSigma_Q\tP-Value\n")
        for snp in snplist:
            fmtstring = "{:s}\t{:d}\t{:d}\t{:g}\t{:g}\t{:g}\t{:g}\t{:g}\n"
            mfile.write(fmtstring.format(snp.rsid, snp.chrom, snp.pos, snp.maf, snp.q, snp.mu, snp.sigma, snp.p))

def prune_region(region, myldict):
    start = time.time()
    sorted_region = sorted(region, key=attrgetter('logp'), reverse=True)
    rejected = collections.defaultdict(lambda: False)
    accepted = []
    for snp in sorted_region:
        if rejected[str(snp.pos)]:
            continue
        accepted.append(snp)
        if myldict[str(snp.pos)]:
            for r in myldict[str(snp.pos)].keys():
                rejected[r] = True
    took = time.time() - start
    return sorted(accepted, key=attrgetter('pos'), reverse=False)

            
def prune_region_GW(region, myldict):
    lonely = list()
    ld_region = collections.defaultdict(list)
    sorted_region = sorted(region, key=attrgetter('logp'), reverse=True)
    rejected = collections.defaultdict(dict)
    for i in range(1,23):
        rejected[i] = collections.defaultdict(lambda: False)
    accepted = []
    accepted_count = collections.defaultdict(int)
    for snp in sorted_region:
        if rejected[snp.chrom][str(snp.pos)]:
            ## add to region
            lead_snp = rejected[snp.chrom][str(snp.pos)]
            ld_region[lead_snp].append(snp.pos)
            # print(snp.chrom, snp.pos, lead_snp)
            continue
        accepted.append(snp)
        if myldict[snp.chrom][str(snp.pos)]:
            for r in myldict[snp.chrom][str(snp.pos)].keys():
                rejected[snp.chrom][r] = snp.rsid
        else:
            #lonely SNP not in LD?
            lonely.append(snp)
    return sort_by_position(accepted), ld_region
    #return sorted(accepted, key=attrgetter('p'), reverse=False), ld_region


def sort_by_position(snps):
    sorted_snps = list()
    for chrom in range(1, 23):
        chrom_snps = [x for x in snps if x.chrom == chrom]
        if len(chrom_snps) > 0:
            sorted_chrom_snps = sorted(chrom_snps, key=attrgetter('pos'), reverse=False)
            sorted_snps += sorted_chrom_snps
    return sorted_snps

def write_ld_region(pruned_snps, ld_regions, outfile):
    with open(outfile, 'w') as outstream:
        for leadsnp in pruned_snps:
            ldr = [str(x) for x in ld_regions[leadsnp.rsid]]
            str_fmt = "{:s}\t{:d}\t{:d}\t{:g}\t{:s}\n".format(leadsnp.rsid, leadsnp.chrom, leadsnp.pos, leadsnp.p, ",".join(ldr))
            outstream.write(str_fmt)
        #for chrm in range(1, 23):
        #    chrm_snps = [x for x in pruned_snps if x.chrom == chrm]
        #    chrm_snps_sorted = sorted(chrm_snps, key=attrgetter('pos'), reverse=False)
        #    for leadsnp in chrm_snps_sorted:

if __name__ == '__main__':

    opts = parse_args()

    print("LDfile:", opts.ldfile)
    LD_gw_dict = dict()
    for chrm in range(1, 23):
        start = time.time()
        print("Loading CHR ", chrm, end="")
        myldict = read_ldfile(opts.ldfile.format(chrm))
        LD_gw_dict[chrm] = myldict
        took = time.time() - start
        print(" - {:g} seconds".format(took))

    pruned_snps = list()
    for i, infile in enumerate(opts.infiles):
        if os.path.exists(infile):
            # prune snps
            print("Prunning ", infile)
            snplist = tejaas_chr(infile)
            pruned_snps, ld_regions = prune_region_GW(snplist, LD_gw_dict)

             # write pruned snps
            #pruned_outfile = infile+".ld_prune"
            pruned_outfile = opts.outfiles[i]
            tejaas_write(pruned_snps, pruned_outfile)
            regions_outfile = os.path.join(os.path.dirname(opts.outfiles[i]), "ld_regions.txt")
            write_ld_region(pruned_snps, ld_regions, regions_outfile)
        else:
            print(infile, "does not exist")
