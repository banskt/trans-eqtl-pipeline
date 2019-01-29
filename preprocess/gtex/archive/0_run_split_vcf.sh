#!/bin/bash

PYTHON="/home/fsimone/myenv/bin/python"
INPUTDIR="/mnt/encrypted/gtex/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/GenotypeFiles/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1"
GTFILE="$INPUTDIR/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz"
TABLEFILE="/mnt/encrypted/gtex/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_VarID_Lookup_Table.txt"

OUTDIR="/home/fsimone/datasets/gtex/genotypes/vcfs_allsamples"
mkdir -p $OUTDIR

OUTPREFIX="$OUTDIR/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_imput_info04_PASS_maf01_HWEp1E6_dbSNP135_ConstrVarIDs_"

$PYTHON split_vcf_in_chr.py --input $GTFILE --rsid-table $TABLEFILE --outprefix $OUTPREFIX --filter-snps
