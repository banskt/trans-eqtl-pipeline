#!/bin/bash
#SBATCH -A cramer
#SBATCH -p em
#SBATCH --qos=short
#SBATCH -t 0-2:00:00
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -J _JOB_NAME
#SBATCH -o _JOB_NAME.out
#SBATCH -e _JOB_NAME.err

PYENV=_PYT_ENV_

export OMP_NUM_THREADS=4

VCFTOOLS=_VCFTOOLS
BGZIP=_BGZIP___
TABIX=_TABIX___
SRCVCF=_SRC_VCF_
ANNOTFILE=_ANNTFILE
OUTFILE=_OUT_FILE
CONVERT_ANNOT_PY=_ANNT_PY_
IMPUTE_MISSING_PY=_IMPT_PY_
VCF_FILTER_PY=_VCFT_PY_
GTVERSION=_GTEX_VER_
CHRM=_CHRM_NUM
MAF_MIN=_MAF_MIN_
REMOVE_INDELS=_RM_INDEL
REMOVE_AMBIGUOUS=_RM_AMBIG
IMPUTE_MISSING=_IMP_MISS

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
if [[ "${GTVERSION}" == "v6p" ]]; then
    sed "18q;d" ${TMPOUT}.recode.vcf | awk -F$'\t' '{for (i = 1; i <= NF; i++) {idmod=gensub(/^(GTEX-[^-]*).*/,"\\1","g", $i); print idmod} }'| awk 'BEGIN { ORS = "\t" } { print }' > ${HEADERFILE}
fi
if [[ "${GTVERSION}" == "v8" ]]; then
    sed "40q;d" ${TMPOUT}.recode.vcf | awk -F$'\t' '{for (i = 1; i <= NF; i++) {idmod=gensub(/^(GTEX-[^-]*).*/,"\\1","g", $i); print idmod} }'| awk 'BEGIN { ORS = "\t" } { print }' > ${HEADERFILE}
fi
${BGZIP} -c ${TMPOUT}.recode.vcf > ${TMPOUT}.recode.vcf.gz
${TABIX} -r ${HEADERFILE} ${TMPOUT}.recode.vcf.gz > ${TMPOUT}.vcf.gz

# remove MAF < 0.1
FILTER_STRING="--maf ${MAF_MIN}"

# remove ambiguous SNPs (optional)
if [ "${REMOVE_AMBIGUOUS}" = "true" ]; then FILTER_STRING="${FILTER_STRING} --remove-ambiguous"; fi
if [[ "${GTVERSION}" == "v8" ]]; then FILTER_STRING="${FILTER_STRING} --maf-metapos 2"; fi

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
