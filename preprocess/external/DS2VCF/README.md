# DS2VCF
Annotate VCF file with dosages estimated from genotype likelihoods


# Installation
To compile `add_dosage` just type `make all`

# Running
The `add_dosage` program assumes a properly formated gzipped vcf file. It only computes dosages
for SNPs and disregards other calls. To compute dosages (DS) it converts phred-scaled likelihoods
(PL), which are assumed to be the final entry for each sample's formatted call. 

    add_dosage [OPTIONS] INPUT.VCF.GZ OUTPUT.VCF.GZ

Right now OPTIONS are disabled.
