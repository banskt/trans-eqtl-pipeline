import argparse
import gzip

def parse_args():

    parser = argparse.ArgumentParser(description='Filter SNPs from a VCF file.')

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

    parser.add_argument('--remove-indels',
                        dest='remove_indels',
                        action='store_true',
                        help='remove indels')

    parser.add_argument('--remove-ambiguous',
                        dest='remove_ambiguous',
                        action='store_true',
                        help='remove ambiguous SNPs')

    parser.add_argument('--qc-pass',
                        dest='vcf_pass_filter',
                        action='store_true',
                        help='remove SNPs without default PASS filter')

    parser.add_argument('--maf',
                        type=float,
                        dest='maf_min',
                        help='remove SNPs with MAF below this value')


    opts = parser.parse_args()
    return opts

if __name__ == '__main__':

    opts = parse_args()

    infile = opts.infile
    outfile = opts.outfile

    remove_indels = opts.remove_indels
    remove_ambiguous = opts.remove_ambiguous
    vcf_pass_filter = opts.vcf_pass_filter
    maf_filter = False
    if opts.maf_min is not None:
        maf_min = opts.maf_min
        maf_filter = True

    SNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

    if remove_ambiguous or remove_indels or vcf_pass_filter or maf_filter:
        nsnps = 0
        nindel = 0
        nambig = 0
        npass = 0
        nmaf = 0
        with gzip.open(infile, 'r') as vcfstream, gzip.open(outfile, 'w') as outstream:
            for line in vcfstream:
                linestrip = line.decode().strip()
                if linestrip[0] == '#':
                    outstream.write(line)
                else:
                    linesplit = linestrip.split("\t")

                    # start filtering
                    deletesnp = False

                    if remove_ambiguous or remove_indels:
                        ref = linesplit[3]
                        alt = linesplit[4]
                        if len(ref) > 1 or len(alt) > 1:
                            if remove_indels:
                                deletesnp = True
                                nindel += 1
                        else:
                            # Ambiguous SNPs are possible only for non-indels
                            if remove_ambiguous:
                                if SNP_COMPLEMENT[ref] == alt:
                                    deletesnp = True
                                    nambig += 1

                    if not deletesnp and vcf_pass_filter:
                        QC = linesplit[6]
                        if QC != "PASS":
                            deletesnp = True
                            npass += 1

                    if not deletesnp and maf_filter:
                        QC_meta = linesplit[7].split(";")
                        maf = float(QC_meta[0].split("=")[1])
                        if maf < maf_min or maf > (1 - maf_min):
                            deletesnp = True
                            nmaf += 1

                    # print SNP if not filtered
                    if not deletesnp:
                        outstream.write(line)
                        nsnps += 1

        if remove_indels:
            print (f"{nindel:d} indels removed")
        if remove_ambiguous:
            print (f"{nambig:d} ambiguous SNPs removed")
        if vcf_pass_filter:
            print (f"{npass:d} SNPs did not pass VCF QC")
        if maf_filter:
            print (f"{nmaf:d} SNPs did not pass MAF criteria")

    else:
        print("Nothing to do. Did you forget specifying the filter?")
