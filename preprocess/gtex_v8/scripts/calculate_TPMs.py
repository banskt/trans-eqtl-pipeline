import numpy as np
import pandas as pd
import gzip
import argparse

def parse_args():

    parser = argparse.ArgumentParser(description='Convert phASER collapsed counts to TPMs')

    parser.add_argument('--in',
                        type=str,
                        dest='infile',
                        metavar='FILE',
                        help='Input read counts')

    parser.add_argument('--out',
                        type=str,
                        dest='outfile',
                        metavar='FILE',
                        help='Output TPMs')

    parser.add_argument('--gtf',
                        type=str,
                        dest='gtffile',
                        metavar='FILE',
                        help='GENCODE GTF file')


    opts = parser.parse_args()
    return opts


opts = parse_args()

### CHECK https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model
### Important for preprocessing the annotation GTF

infile="/cbscratch/franco/datasets/gtex_v8/gencode.v26.collapsed.genes.gtf"
outfile="/cbscratch/franco/datasets/gtex_v8/gencode.v26.coding_lengths"

l = 0
genedict = collections.defaultdict(dict)
lead_transcripts = dict()
with open(infile) as instream:
    for line in instream:
        if line.startswith("#"):
            continue
        arr = line.strip().split("\t")
        chrm = arr[0]
        etype = arr[2]
        start = int(arr[3])
        end = int(arr[4])
        annots = arr[8].split(";")
        gene_id = annots[0].split()[1][1:-1]
        if etype == "gene":
            genedict[gene_id]["length"] = end - start + 1
        if etype == "exon":
            for aa in annots:
                if aa.startswith(" transcript_id"):
                    t_id = aa.split()[1][1:-1]
                    if t_id in genedict[gene_id]:
                        genedict[gene_id][t_id] += end - start + 1
                    else:
                        genedict[gene_id][t_id] = end - start + 1
                        lead_transcripts[gene_id] = t_id
                    break
        l += 1
#         if l > 20:
#             break

with open(outfile, 'w') as outstream:
    for i, gene in enumerate(lead_transcripts.keys()):
        t_id = lead_transcripts[gene]
        t_length = genedict[gene][t_id]
        outstream.write("{:s}\t{:d}\n".format(gene, t_length))

### Read gene lengths file for gencode v26
lengthsfile = outfile="/cbscratch/franco/datasets/gtex_v8/gencode.v26.coding_lengths"
gene_len_dict = dict()

with open(lengthsfile) as instream:
    for line in instream:
        arr = line.rstrip().split()
        gene_len_dict[arr[0].split(".")[0]] = int(arr[1])

ix = list()
genenames = [x.split(".")[0] for x in list(reads.Name)]
for name in genenames:
    if name in gene_len_dict:
        ix.append(True)
    else:
        ix.append(False)

geneix = np.where(ix)[0]
filtered_genenames = [genenames[i] for i in geneix]
genelengths = np.array([gene_len_dict[x] for x in filtered_genenames])

##############################
# Read count file for phASER expression matrix
##############################

gtex_v8_reads = opts.infile #"/cbscratch/franco/datasets/gtex_v8/expression/phASER_GTEx_v8_matrix_collapsed_counts_new.txt.gz"
outfile = opts.outfile # "/cbscratch/franco/datasets/gtex_v8/expression/TPMs_phASER_GTEx_v8.matrix.txt"
reads_phaser = pd.read_csv(gtex_v8_reads, header=0, sep="\t")





##############################
# Calculate TPM values for phASER expression matrix
##############################

filtered_reads = reads_phaser.iloc[np.array(ix)]
TPM_matrix = np.zeros((filtered_reads.shape[0], filtered_reads.shape[1]-1))
for i in range(1, filtered_reads.shape[1]):
    rpk = filtered_reads.iloc[:, i] * 1000 / genelengths
    sf = np.sum(rpk) / 1000000
    tpm = rpk / sf
    TPM_matrix[:,i-1] = tpm
    if i%1000 == 0:
        print("Processed {:d} samples".format(i))

tpm_df = pd.DataFrame(TPM_matrix, index=filtered_reads.name, columns=filtered_reads.columns[1:])
tpm_df.to_csv(outfile, header=True, sep="\t")