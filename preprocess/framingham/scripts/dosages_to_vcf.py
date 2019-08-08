import gzip
import argparse
import collections
import re
import numpy as np

def parse_args():

    parser = argparse.ArgumentParser(description='Convert dosages to vcf.gz file')

    parser.add_argument('--dosagein',
                        type=str,
                        dest='in_file',
                        metavar='FILE',
                        help='input dosage file of a single chromosome')

    parser.add_argument('--outfile',
                        type=str,
                        dest='out_file',
                        metavar='FILE',
                        help='output prefix for chromosome file, including output directory')

    parser.add_argument('--samplefile',
                        type=str,
                        dest='sample_file',
                        metavar='STR',
                        help='sample file in oxford format')

    parser.add_argument('--notdosage',
                        dest='not_dosage',
                        action='store_true', 
                        help='instead of dosage (0-2) is a combination of 3N fields for each sample (for Cardiogenics)')

    parser.set_defaults(not_dosage=False)

    parser.add_argument('--callthres',
                        dest='call_thres',
                        type=float, 
                        help='GT calling threshold for dosages',
                        default=0.1)

    opts = parser.parse_args()
    return opts

def read_samples(samplefile):
    with open(samplefile, 'r') as samfile:
        sample = 0
        samplenames = list()
        next(samfile)
        next(samfile)
        for line in samfile:
            if re.search('^#', line):
                continue
            sample += 1
            samplenames.append(line.strip().split()[0])
    nsample = sample
    return nsample, samplenames

if __name__ == '__main__':
    
    args = parse_args()

    filepath = args.in_file
    outfilepath = args.out_file
    samplefilepath = args.sample_file

    nsamples, samplenames = read_samples(samplefilepath)

    nloci = 0
    linenum = 0

    if args.not_dosage:
        meta_columns = 5
    else:
        meta_columns = 6

    vcffile = gzip.open(outfilepath, 'w')

    vcffile.write("##fileformat=VCFv4.1\n".encode('utf-8'))
    vcffile.write("##fileDate=201907\n".encode('utf-8'))
    vcffile.write("##source=Converted from dosages\n".encode('utf-8'))
    vcffile.write("##reference=Homo_sapiens_assembly19_v1\n".encode('utf-8'))
    vcffile.write("##INFO=<ID=EXP_FREQ_A1,Number=.,Type=Float,Description=\"Expected frequency of allele coded 1\">\n".encode('utf-8'))
    vcffile.write("##FORMAT=<ID=GT:,Number=1,Type=String,Description=\"Best Guessed Genotype with posterior probability threshold of 0.9\">\n".encode('utf-8'))
    vcffile.write("##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">\n".encode('utf-8'))
    vcffile.write("##FILTER=<ID=maf00,Description=\"allele frequency between 0\% and 100\%\">\n".encode('utf-8'))
    sample_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(samplenames)+"\n"
    vcffile.write(sample_header.encode('utf-8'))

    info_string = "\t.\tPASS\tEXP_FREQ_A1={:f}\tGT:DS\t"

    with gzip.open(filepath, 'r') as filereader:
        for snpline in filereader:
            nloci += 1
            mline = snpline.decode('utf-8').rstrip().split()

            if args.not_dosage:
                ngenotypes = (len(mline) - meta_columns) / 3
            else:
                ngenotypes = len(mline) - meta_columns

            if float(ngenotypes).is_integer():
                if ngenotypes != nsamples:
                    print('Number of samples differ from genotypes')
                    raise;
            else:
                print('Number of columns in genotype frequencies not divisible by 3')
                raise

            if args.not_dosage:
                gt_freqs = np.array([float(x) for x in mline[meta_columns:]])
                indsAA = np.arange(0,nsamples) * 3
                indsAB = indsAA + 1
                indsBB = indsAB + 1
                snp_dosage = 2*gt_freqs[indsBB] + gt_freqs[indsAB] # [AA, AB, BB] := [0, 1, 2]
            else:
                snp_dosage = np.array([float(x) for x in mline[meta_columns:]])

            maf = sum(snp_dosage) / 2 / len(snp_dosage)
            try:                               ######## change to get the chrom numberfrom gtfile
               chrom = int(mline[0])
            except:
               chrom = 0

            bp_pos     = int(mline[2])
            varid      = mline[1]
            ref_allele = mline[3]
            alt_allele = mline[4]
            
            GT_list = list()
            margin = args.call_thres
            for x in snp_dosage:
                if x > (2.0 - margin):
                    GT_list.append("1/1")
                elif x > (1.0 - margin) and x < (1.0 + margin):
                    GT_list.append("0/1")
                elif x < margin:
                    GT_list.append("0/0")
                else:
                    GT_list.append("./.")

            vcffile.write("{:d}\t{:d}\t{:s}\t{:s}\t{:s}".format(chrom, bp_pos, varid, ref_allele, alt_allele).encode('utf-8'))
            vcffile.write(info_string.format(maf).encode('utf-8'))
            dosage_str = "\t".join(["{:s}:{:g}".format(y,x) for y,x in zip(GT_list, snp_dosage)]) + "\n"
            vcffile.write(dosage_str.encode('utf-8'))