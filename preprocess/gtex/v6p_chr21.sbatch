#!/bin/bash
#SBATCH -p medium
#SBATCH --qos=short
#SBATCH -t 0-2:00:00
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -C scratch
#SBATCH -A all
#SBATCH -J _JOB_NAME
#SBATCH -o _JOB_NAME.out
#SBATCH -e _JOB_NAME.err


PYENV=/usr/users/sbanerj/miniconda3/envs/py36/bin/python

export OMP_NUM_THREADS=4

VCFTOOLS=/usr/users/sbanerj/packages/vcftools/vcftools-v0.1.15/bin/vcftools
BGZIP=/usr/users/sbanerj/packages/htslib/htslib-1.4.1/bin/bgzip
TABIX=/usr/users/sbanerj/packages/htslib/htslib-1.4.1/bin/tabix
SRCVCF=/cbscratch/sbanerj/data/GTEx/v6p/genotype/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz
ANNOTFILE=/cbscratch/sbanerj/trans-eqtl/input/gtex_v6p/genotype/all_samples/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs_chr21.annot
OUTFILE=/cbscratch/sbanerj/trans-eqtl/input/gtex_v6p/genotype/all_samples/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs_maf01_noindels_noambig_chr21
CONVERT_ANNOT_PY=/usr/users/sbanerj/trans-eQTL/dev-pipeline/scripts/vcf_change_annot.py
IMPUTE_MISSING_PY=/usr/users/sbanerj/trans-eQTL/dev-pipeline/scripts/vcf_impute_missingGT.py
VCF_FILTER_PY=/usr/users/sbanerj/trans-eQTL/dev-pipeline/scripts/vcf_filter.py
CHRM=21
MAF_MIN=0.01
REMOVE_INDELS=true
REMOVE_AMBIGUOUS=true
IMPUTE_MISSING=true

MAF_MAX=$( echo $MAF_MIN | awk '{m = 1 - $1; printf "%g", m}' )
OUTDIR="$(dirname "${OUTFILE}")"
RAND=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
TMPOUT="${OUTDIR}/tmp_${RAND}_chr${CHRM}"
HEADERFILE="${OUTDIR}/tmp_${RAND}_chr${CHRM}_newheader.txt"


# split and GTEx filter and remove indels (optional)
VCF_FILTER_OPTIONS="--remove-filtered-all"
if [ "${REMOVE_INDELS}" = "true" ]; then
    VCF_FILTER_OPTIONS="${VCF_FILTER_OPTIONS} --remove-indels"
fi
${VCFTOOLS} --gzvcf ${SRCVCF} --chr ${CHRM} ${VCF_FILTER_OPTIONS} --recode --recode-INFO-all --out ${TMPOUT}

# change header
sed "18q;d" ${TMPOUT}.recode.vcf | awk -F$'\t' '{for (i = 1; i <= NF; i++) {idmod=gensub(/^(GTEX-[^-]*).*/,"\\1","g", $i); print idmod} }'| awk 'BEGIN { ORS = "\t" } { print }' > ${HEADERFILE}
${BGZIP} -c ${TMPOUT}.recode.vcf > ${TMPOUT}.recode.vcf.gz
${TABIX} -r ${HEADERFILE} ${TMPOUT}.recode.vcf.gz > ${TMPOUT}.vcf.gz

# remove MAF < 0.1
FILTER_STRING="--maf ${MAF_MIN}"

# remove ambiguous SNPs (optional)
if [ "${REMOVE_AMBIGUOUS}" = "true" ]; then
    FILTER_STRING="${FILTER_STRING} --remove-ambiguous"
fi
${PYENV} ${VCF_FILTER_PY} --input ${TMPOUT}.vcf.gz --out ${TMPOUT}_pyfilter.vcf.gz ${FILTER_STRING}

# convert annotation
${PYENV} ${CONVERT_ANNOT_PY} --input ${TMPOUT}_pyfilter.vcf.gz --out ${OUTFILE}.vcf.gz --annot ${ANNOTFILE}
gunzip ${OUTFILE}.vcf.gz
${BGZIP} ${OUTFILE}.vcf # ${PYENV} creates gunzipped files, not bgzipped, Tabix needs bgzipped.
${TABIX} -f ${OUTFILE}.vcf.gz

# impute missing genotype (optional)
if [ "${IMPUTE_MISSING}" = "true" ]; then
    ${PYENV} ${IMPUTE_MISSING_PY} --input ${OUTFILE}.vcf.gz --out ${OUTFILE}_WARNimputedmissing.vcf.gz
    gunzip ${OUTFILE}_WARNimputedmissing.vcf.gz
    ${BGZIP} ${OUTFILE}_WARNimputedmissing.vcf
fi

# remove temp files
rm -rf ${TMPOUT}* ${HEADERFILE}
