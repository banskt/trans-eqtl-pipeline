import argparse
import gzip
import numpy as np
import pandas as pd
import os 

def parse_args():

    parser = argparse.ArgumentParser(description='Filter significant trans-eQTLs from RR output.')

    parser.add_argument('--indir',
                        type=str,
                        dest='indir',
                        metavar='FILE',
                        help='directory of chromosomes with rr.txt and gene_snp_list.txt compiled chunks file')

    parser.add_argument('--c',
                        type=float,
                        dest='cutoff',
                        default=5e-8,
                        help='significance cutoff for trans-eqtls')

    parser.add_argument('--outdir',
                        type=str,
                        dest='outdir',
                        metavar='FILE',
                        help='output directory')    

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':

    opts = parse_args()

    chrms = [ "chr"+str(x) for x in np.arange(1,23)]
    with open(opts.outdir+"/trans_eqtls_{:s}.txt".format(str(opts.cutoff)), 'w') as outstream:
        with open(opts.outdir+"/target_genes_{:s}.txt".format(str(opts.cutoff)), 'w') as outstream2:
            with open(opts.outdir+"/target_genes_knn_{:s}.txt".format(str(opts.cutoff)), 'w') as outstream_knn:
                with open(opts.outdir+"/snps_list.txt", 'w') as outstream3:
                    with open(opts.outdir+"/target_genes_FDR_{:s}.txt".format(str(opts.cutoff)), 'w') as outstream_fdr:
                        for chrm in chrms:
                            rr_infile = os.path.join(opts.indir, chrm, "rr.txt")
                            gene_infile = os.path.join(opts.indir, chrm, "gene_snp_list.txt")
                            knn_gene_infile = os.path.join(opts.indir, chrm, "gene_snp_list_knn.txt")
                            fdr_gene_infile = os.path.join(opts.indir, chrm, "gene_snp_list_FDR.txt")
                            if os.path.exists(rr_infile) and os.path.exists(gene_infile):
                                df = pd.read_csv(rr_infile, header=0, sep="\t")
                                df2 = pd.read_csv(gene_infile, header=0, sep="\t")
                                snps = df[ ["ID", "MAF"] ]
                                outdf = df[ df["P"] <= opts.cutoff ]                  
                                outdf = outdf[["ID", "CHR", "Pos", "MAF", "Q", "Mu", "Sigma", "P"]]
                                out_targets = df2[ (df2['snpid'].isin(outdf["ID"])) & (df2['pval']<0.01) ]

                                if os.path.exists(knn_gene_infile):
                                    df3 = pd.read_csv(knn_gene_infile, header=0, sep="\t")
                                    knn_out_targets = df3[ (df3['snpid'].isin(outdf["ID"])) & (df3['pval']<0.01) ]

                                if os.path.exists(fdr_gene_infile):
                                    df4 = pd.read_csv(fdr_gene_infile, header=0, sep="\t")
                                    fdr_out_targets = df4[ (df4['snpid'].isin(outdf["ID"])) & (df4['pval']<0.01) ]

                                if chrm == "chr1":
                                    outdf.to_csv(outstream, sep="\t", header=True, index=False, float_format='%g')
                                    out_targets.to_csv(outstream2, sep="\t", header=True, index=False, float_format='%g')
                                    if os.path.exists(knn_gene_infile):
                                        knn_out_targets.to_csv(outstream_knn, sep="\t", header=True, index=False, float_format='%g')
                                    if os.path.exists(fdr_gene_infile):
                                        fdr_out_targets.to_csv(outstream_fdr, sep="\t", header=True, index=False, float_format='%g')
                                    snps.to_csv(outstream3, sep="\t", header=True, index=False, float_format='%g')
                                else:
                                    outdf.to_csv(outstream, sep="\t", header=False, index=False, float_format='%g')
                                    out_targets.to_csv(outstream2, sep="\t", header=False, index=False, float_format='%g')
                                    if os.path.exists(knn_gene_infile):
                                        knn_out_targets.to_csv(outstream_knn, sep="\t", header=False, index=False, float_format='%g')
                                    if os.path.exists(fdr_gene_infile):
                                        fdr_out_targets.to_csv(outstream_fdr, sep="\t", header=False, index=False, float_format='%g')
                                    snps.to_csv(outstream3, sep="\t", header=False, index=False, float_format='%g')
                            else:
                                print("rr.txt or gene_snp_list.txt does not exist ->", rr_infile)
                                raise