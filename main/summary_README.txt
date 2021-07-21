# This directory contains the results of TEJAAS on GTEx version 8.
# Each folder corresponds to one of the 49 tissues available. The abbreviations definitions can be found at the end of this file.

# Each folder contains the following files:

- trans_eqtls.txt
Contains the list of all trans-eQTLs SNPs which obtained a RR score p-value < 5e-8. 
It has 8 columns:
  - ID: SNP variant id. It can be converted to rsid using the lookup table grom the GTEx Portal (https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz)
  - CHR: chromosome number
  - Pos: genomic position in that chromosome
  - MAF: Minor allele frequency for that SNP in the set of samples which have gene expression available in that tissue. The same SNP might have different MAF across tissues, since not all samples have gene expression available for all tissues.
  - Q: Reverse Regression score (RR-score), also called Q-rev in our manuscript. 
  - Mu: analytically estimated mean of the RR-score distribution
  - Sigma: analytically estimated standard deviation of the RR-score distribution
  - P: p-value for Q, calculated using the CDF of a Normal distribution N(Mu, Sigma)

- trans_eqtls_ldpruned.txt
List of all the lead trans-eQTLs SNPs from the file trans_eqtls.txt. LD pruned was performed a window size of 200Kb and discarded all SNPs with R^2 > 0.5

- ld_regions.txt
Contains a list of all the LD regiones with significant trans-eQTLs discovered with TEJAAS.
It has 5 columns:
  - Lead SNP variant id.
  - Chromosome number
  - Genomic position of lead SNP
  - P-value for RR-score of lead SNP
  - List of genomic positions for significant SNPs in LD with the lead SNP.

- target_genes.txt
List of all possible target genes for each SNP in the trans_eqtls.txt file. Each row is a SNP-gene pair
It has 4 columns:
  - geneid: Ensembl gene id
  - snpid: SNP variant id
  - pval: P-value of the simple linear regression association of that SNP-gene pair
  - adj_pval: Adjusted p-value for Benjamini Hochberg (BH) procedure with an FDR of 50%. BH was calculated for each SNP and all its possible target genes. FDR values are not comparable across SNPs.


# List of Abbreviations
as      Adipose - Subcutaneous
av      Adipose - Visceral (Omentum)
ag      Adrenal Gland
aa      Artery - Aorta
ac      Artery - Coronary
at      Artery - Tibial
bam     Brain - Amygdala
ban     Brain - Anterior cingulate cortex (BA24)
bca     Brain - Caudate (basal ganglia)
bceh    Brain - Cerebellar Hemisphere
bce     Brain - Cerebellum
bco     Brain - Cortex
bfr     Brain - Frontal Cortex (BA9)
bhi     Brain - Hippocampus
bhy     Brain - Hypothalamus
bnu     Brain - Nucleus accumbens (basal ganglia)
bpu     Brain - Putamen (basal ganglia)
bsp     Brain - Spinal cord (cervical c-1)
bsu     Brain - Substantia nigra
br      Breast - Mammary Tissue
ebv     Cells - EBV-transformed lymphocytes
fib     Cells - Cultured fibroblasts
cols    Colon - Sigmoid
colt    Colon - Transverse
esog    Esophagus - Gastroesophageal Junction
esom    Esophagus - Mucosa
esomu   Esophagus - Muscularis
haa     Heart - Atrial Appendage
hlv     Heart - Left Ventricle
kc      Kidney - Cortex
liv     Liver
lu      Lung
msg     Minor Salivary Gland
ms      Muscle - Skeletal
nt      Nerve - Tibial
pan     Pancreas
pit     Pituitary
snse    Skin - Not Sun Exposed (Suprapubic)
sse     Skin - Sun Exposed (Lower leg)
si      Small Intestine - Terminal Ileum
spl     Spleen
sto     Stomach
thy     Thyroid
wb      Whole Blood
ov      Ovary
pro     Prostate
tes     Testis
ut      Uterus
va      Vagina
