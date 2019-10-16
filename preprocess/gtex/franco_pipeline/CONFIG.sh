#!/bin/bash

TISSUEFILE="tissues.table.txt"
ENV="$HOME/opt/miniconda/3/envs/env3.6/bin/python"

BASEDIR="/cbscratch/franco/datasets/gtex"

# REQUIRED FILES!! download them from GTEx (for pheno file you need access)
RPKMFILE="GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz"
SAMPLE_PHENOFILE="phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt"
SUBJECT_PHENOFILE="phs000424.v6.pht002742.v6.p1.c1.GTEx_Subject_Phenotypes.GRU.txt"
READSFILE="GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz"
DONORFILE="$BASEDIR/donor_ids.fam"


RPKMOUTDIR="$BASEDIR/expression/rpkms"              # tissue specific rpkms output dir
EXPROUTDIR="$BASEDIR/expression/normalized_expr"    # normalized, uncorrected expr dir

GTEXCOVDIR="$BASEDIR/expression/GTEx_Analysis_v6p_eQTL_covariates"  # dir for covariates downloaded from GTEx (public access)
COVDIR="$BASEDIR/expression/PEER_outputs"                        # PEER covariates outdir
RUNPEER=false
MAXNCOV=35                     # max number of Nº of hidden confounders for PEER

LMOUTDIR="$BASEDIR/expression/norm_lmcorrected" # outdir of linear model corrected expressions (PC+platform+sex, no peer!)

# If you want to correct for PEER, use _residuals output directly from step 2 script.

VCFDIR="$BASEDIR/genotypes/vcfs_allsamples"
PLINK2="$HOME/bin/plink2"
