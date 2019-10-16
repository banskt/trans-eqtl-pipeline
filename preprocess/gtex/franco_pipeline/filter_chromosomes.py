import sys
sys.path.append("../")
from iotools.readOxford import ReadOxford
import argparse

def parse_args():

    parser = argparse.ArgumentParser(description='Chromosome filtering')
    parser.add_argument('--oxf',
                        type=str,
                        dest='oxf',
                        metavar='DIR',
                        help='input genotype file')
    parser.add_argument('--fam',
                        type=str,
                        dest='fam',
                        metavar='DIR',
                        help='Samples input file')

    parser.add_argument('--chrom',
                        type=int,
                        dest='chrom',
                        metavar='CHROM',
                        help='Chromosome number to use')

    parser.add_argument('--out',
                        type=str,
                        dest='out',
                        metavar='OUTFILE',
                        help='Outfile to write filtered genotype')
    parser.add_argument('--dset',
                        type=str,
                        dest='dataset',
                        metavar='DATASET',
                        help='gtex or cardiogenics')

    opts = parser.parse_args()
    return opts

opts = parse_args()

print("Chr {:d}".format(opts.chrom))
oxf_file = opts.oxf # "/cbscratch/franco/datasets/cardiogenics/genotypes/CG_{:d}.imputed.gz".format(chrom)
outfile  = opts.out # "/cbscratch/franco/datasets/cardiogenics/genotypes/prefiltered/CG_dosages_filtered_{:d}.imputed.gz".format(chrom)
fam_file = opts.fam # "/cbscratch/franco/datasets/cardiogenics/genotypes/CG.sample"

if opts.dataset == "cardiogenics":
    oxf = ReadOxford(oxf_file, fam_file, isdosage=False, chrom=opts.chrom) #should be self.args.isdosage and self.args.oxf_columns
if opts.dataset == "gtex":
    oxf = ReadOxford(oxf_file, fam_file, isdosage=True, chrom=opts.chrom) #should be self.args.isdosage and self.args.oxf_columns
dosage = oxf.dosage
gt_donor_ids = oxf.samplenames
snpinfo = oxf.snpinfo

oxf.write_dosages(outfile, format="gtex", filter_gt=True)
