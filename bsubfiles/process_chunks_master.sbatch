#!/bin/bash
#SBATCH -p medium
#SBATCH --qos=short
#SBATCH -t 0-2:00:00
#SBATCH -n 1
#SBATCH -J _JOB_NAME
#SBATCH -o _JOB_NAME.out
#SBATCH -e _JOB_NAME.err

TBASE=__TBASE__
TSHORT=__TSHORT__

source /usr/users/sbanerj/trans-eQTL/dev-pipeline/main/utils/tejaas_chunk_reduce.new

CHRNUMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
OUTDIR="/usr/users/sbanerj/trans_eqtl_results/gtex_v8"
EXPRESSIONFILE="/usr/users/sbanerj/gtex_v8/expression/gtexv8_${TSHORT}_raw_std_protein_coding.txt"
METHOD_VARIANT="permnull_sb0.1_knn"
TEJAAS_PREPROC_STR="raw_std"

RANDSTRING=$( cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1 )

## chromosomes in separate outer loop to prevent reading VCF file multiple times.
NCHUNK_IN_CHRM=(0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0) # initialize with 22 elements
for CHRM in ${CHRNUMS}; do 
    CHRINDX=$((CHRM - 1))
    SPECIFIC_OUTDIR="/usr/users/sbanerj/trans-eQTL/dev-pipeline/jobsubs/gtex_v8/as/tejaas/raw_std/permnull_sb0.1_knn/chr${CHRM}"
    NCHUNK=$( ls -l ${SPECIFIC_OUTDIR}/*.sbatch | wc -l )
    NCHUNK_IN_CHRM[${CHRINDX}]=${NCHUNK}
done


OUTDIR_DATA="${OUTDIR}/${TSHORT}"

if [ ! -z "${EXPRESSIONFILE}" ]; then

    for CHRM in ${CHRNUMS}; do
        CHRINDX=$((CHRM - 1))
        NCHUNK=${NCHUNK_IN_CHRM[${CHRINDX}]}
        echo "chr${CHRM} --> ${NCHUNK} chunks"
        echo "============================"
        tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas/${TEJAAS_PREPROC_STR}/${METHOD_VARIANT}/chr${CHRM}" ${NCHUNK}
    done

    ## collect the trans-eQTLs
    __COLRES="/usr/users/sbanerj/trans_eqtl_results/gtex_v8_tejaas_${METHOD_VARIANT}/${TSHORT}"
    if [ ! -d ${__COLRES} ]; then mkdir -p ${__COLRES}; fi

    __RROUT="${OUTDIR_DATA}/tejaas/${TEJAAS_PREPROC_STR}/${METHOD_VARIANT}/trans_eqtls.txt"
    __TGOUT="${OUTDIR_DATA}/tejaas/${TEJAAS_PREPROC_STR}/${METHOD_VARIANT}/target_genes.txt"
    echo -e "ID\tCHROM\tPOS\tP-VALUE" > ${__RROUT}
    echo -e "ID\tTARGET_GENE\tP-VALUE" > ${__TGOUT}
    echo -ne "Creating list of trans-eQTLs and target genes from "
    for CHRM in ${CHRNUMS}; do
        echo -ne "${CHRM} "
        __TMPDIR="${OUTDIR_DATA}/tejaas/${TEJAAS_PREPROC_STR}/${METHOD_VARIANT}/chr${CHRM}/tmp_${RANDSTRING}"
        if [ -d ${__TMPDIR} ]; then rm -rf ${__TMPDIR}; fi; mkdir -p ${__TMPDIR};
        __CHROUTDIR="${OUTDIR_DATA}/tejaas/${TEJAAS_PREPROC_STR}/${METHOD_VARIANT}/chr${CHRM}"
        cat ${__CHROUTDIR}/rr.txt | awk '$6 < 5e-8 {print}' | awk -v CHRMNUM="${CHRM}" '{printf "%s\t%d\t%d\t%g\n",  $1, CHRMNUM, $2, $6}' > ${__TMPDIR}/${TSHORT}_chr${CHRM}.txt
        cut -f1 ${__TMPDIR}/${TSHORT}_chr${CHRM}.txt > ${__TMPDIR}/${TSHORT}_chr${CHRM}_teqtl.txt
        cat ${__TMPDIR}/${TSHORT}_chr${CHRM}.txt >> ${__RROUT}
        grep -Fwf ${__TMPDIR}/${TSHORT}_chr${CHRM}_teqtl.txt ${__CHROUTDIR}/gene_snp_list.txt | awk '{printf "%s\t%s\t%g\n", $2, $1, $3}' >> ${__TGOUT}
        rm -rf ${__TMPDIR}
    done
    echo "Done."

    cp ${__RROUT} ${__COLRES}/
    cp ${__TGOUT} ${__COLRES}/

fi
