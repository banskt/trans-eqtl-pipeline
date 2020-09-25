import numpy as np
import pandas as pd
import gzip
import argparse
import collections
import re 

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

infile=opts.gtffile 
outfile="/cbscratch/franco/datasets/geuvadis/gencode.v12.geuvadis.coding_lengths"

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

with open(outfile, 'w') as outstream:
    for i, gene in enumerate(lead_transcripts.keys()):
        t_id = lead_transcripts[gene]
        t_length = genedict[gene][t_id]
        outstream.write("{:s}\t{:d}\n".format(gene, t_length))

### Read gene lengths file for gencode v26
lengthsfile = outfile  # this is the same file from above
gene_len_dict = dict()

with open(lengthsfile) as instream:
    for line in instream:
        arr = line.rstrip().split()
        gene_len_dict[arr[0].split(".")[0]] = int(arr[1])

# ix = list()
# genenames = [x.split(".")[0] for x in list(reads.Name)]
# for name in genenames:
#     if name in gene_len_dict:
#         ix.append(True)
#     else:
#         ix.append(False)

# geneix = np.where(ix)[0]
# filtered_genenames = [genenames[i] for i in geneix]


##############################
# Read count file of expression matrix
##############################

geu_reads = opts.infile # Counts file
outfile = opts.outfile  

reads_phaser = pd.read_csv(geu_reads, header=0, sep="\t")
reads_phaser.columns = [x.strip() for x in reads_phaser.columns]
genenames = list(reads_phaser.iloc[:,0])
genelengths = np.array([gene_len_dict[x.split(".")[0]] for x in genenames])

##############################
# Calculate TPM values for phASER expression matrix
##############################

TPM_matrix = np.zeros((reads_phaser.shape[0], reads_phaser.shape[1]-1))
for i in range(1, reads_phaser.shape[1]):
    rpk = reads_phaser.iloc[:, i] * 1000 / genelengths
    sf = np.sum(rpk) / 1000000
    tpm = rpk / sf
    TPM_matrix[:,i-1] = tpm
    if i%1000 == 0:
        print("Processed {:d} samples".format(i))

tpm_df = pd.DataFrame(TPM_matrix, index=reads_phaser.iloc[:,0], columns=reads_phaser.columns[1:])
tpm_df.to_csv(outfile, header=True, sep="\t")