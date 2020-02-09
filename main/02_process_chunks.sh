#!/bin/bash

if [ "$PS1" ]; then echo -e "This script cannot be sourced. Use \"${BASH_SOURCE[0]}\" instead." ; return ; fi

CONFIGFILE=$1

if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ./02_process_chunks.sh CONFIGFILE"
    exit 1
fi

source ${CONFIGFILE}
source PATHS
source ${UTILSDIR}/unset_vars
source ${UTILSDIR}/tejaas_chunk_reduce.new

source ${DATALOAD}

RANDSTRING=$( cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1 )

## chromosomes in separate outer loop to prevent reading VCF file multiple times.
NCHUNK_IN_CHRM=(0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0) # initialize with 22 elements
for CHRM in ${CHRNUMS}; do 
    CHRINDX=$((CHRM - 1))
    SPECIFIC_OUTDIR="${JOBSUBDIR}/ag/tejaas/raw_std/permnull_sb0.002_knn/chr${CHRM}"
    NCHUNK=$( ls -l ${SPECIFIC_OUTDIR}/*.sbatch | wc -l )
    #NCHUNK=$( expected_nchunk ${GENO_FMT} ${CHRM} ${MAX_NSNP_PERJOB} )
    NCHUNK_IN_CHRM[${CHRINDX}]=${NCHUNK}
    echo "Chr${CHRM} was split into ${NCHUNK} jobs."
done

while IFS='' read -r LINE || [ -n "$LINE" ]; do
    if [[ $LINE =~ ^[^\#] ]]; then

        # Get the name of the tissue
        TFULL=$( echo "${LINE}" | cut -f 1 )
        TSHORT=$( echo "${LINE}" | cut -f 2 )
        TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

        OUTDIR_DATA="${OUTDIR}/${TSHORT}"
        GX_TISSUE_FMT=${EXPR_FMT/\[TISSUE\]/${TSHORT}}
        EXPRESSIONFILE=${GX_TISSUE_FMT/\[PREPROC_STRING\]/${TEJAAS_PREPROC_STR}}

        echo ${OUTDIR_DATA}

        if [ ! -z "$EXPRESSIONFILE" ]; then

            echo "Processing chunks for ${TSHORT} (${TBASE})"

            for CHRM in ${CHRNUMS}; do
                CHRINDX=$((CHRM - 1))
                NCHUNK=${NCHUNK_IN_CHRM[${CHRINDX}]}
                echo "chr${CHRM} --> ${NCHUNK} chunks"
                echo "============================"
                for NULL in ${TEJAAS_NULL}; do
                    if [ ${NULL} == "perm" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_PERM}; fi
                    if [ ${NULL} == "maf" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_MAF}; fi
                    for SBETA in ${TEJAAS_SIGMA_BETA}; do
                        METHOD_VARIANT="${NULL}null_sb${SBETA}"
                        if [ "${TEJAAS_KNN}" = "true" ]; then METHOD_VARIANT="${METHOD_VARIANT}_knn"; fi
                        if [ "${bTejaas}" = "true" ];    then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas/${TEJAAS_PREPROC_STR}/${METHOD_VARIANT}/chr${CHRM}" $NCHUNK; fi
                        if [ "${bTjsRandom}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_DATA}/tejaas_rand/${TEJAAS_PREPROC_STR}/${METHOD_VARIANT}/chr${CHRM}" $NCHUNK; fi
                    done
                done
            done

            ## collect the trans-eQTLs
            #if [ "${bTejaas}" = "true" ]; then
            #    METHOD_VARIANT="${TEJAAS_NULL}null_sb${TEJAAS_SIGMA_BETA_PERM}_knn"

            #    __COLRES="/usr/users/sbanerj/trans_eqtl_results/gtex_v8_tejaas_${METHOD_VARIANT}/${TSHORT}"
            #    if [ ! -d ${__COLRES} ]; then mkdir -p ${__COLRES}; fi

            #    __RROUT="${OUTDIR_DATA}/tejaas/${TEJAAS_PREPROC_STR}/${METHOD_VARIANT}/trans_eqtls.txt"
            #    __TGOUT="${OUTDIR_DATA}/tejaas/${TEJAAS_PREPROC_STR}/${METHOD_VARIANT}/target_genes.txt"
            #    echo -e "ID\tCHROM\tPOS\tP-VALUE" > ${__RROUT}
            #    echo -e "ID\tTARGET_GENE\tP-VALUE" > ${__TGOUT}
            #    echo -ne "Creating list of trans-eQTLs and target genes from "
            #    for CHRM in ${CHRNUMS}; do
            #        echo -ne "${CHRM} "
            #        __TMPDIR="${OUTDIR_DATA}/tejaas/${TEJAAS_PREPROC_STR}/${METHOD_VARIANT}/chr${CHRM}/tmp_${RANDSTRING}"
            #        if [ -d ${__TMPDIR} ]; then rm -rf ${__TMPDIR}; fi; mkdir -p ${__TMPDIR};
            #        __CHROUTDIR="${OUTDIR_DATA}/tejaas/${TEJAAS_PREPROC_STR}/${METHOD_VARIANT}/chr${CHRM}"
            #        cat ${__CHROUTDIR}/rr.txt | awk '$6 < 5e-8 {print}' | awk -v CHRMNUM="${CHRM}" '{printf "%s\t%d\t%d\t%g\n",  $1, CHRMNUM, $2, $6}' > ${__TMPDIR}/${TSHORT}_chr${CHRM}.txt
            #        cut -f1 ${__TMPDIR}/${TSHORT}_chr${CHRM}.txt > ${__TMPDIR}/${TSHORT}_chr${CHRM}_teqtl.txt
            #        cat ${__TMPDIR}/${TSHORT}_chr${CHRM}.txt >> ${__RROUT}
            #        grep -Fwf ${__TMPDIR}/${TSHORT}_chr${CHRM}_teqtl.txt ${__CHROUTDIR}/gene_snp_list.txt | awk '{printf "%s\t%s\t%g\n", $2, $1, $3}' >> ${__TGOUT}
            #        rm -rf ${__TMPDIR}
            #    done
            #    echo "Done."

            #    cp ${__RROUT} ${__COLRES}/
            #    cp ${__TGOUT} ${__COLRES}/
            #fi

        fi
    fi
done < ${TISSUEFILE}

unset_vars ${CONFIGFILE}
unset_vars PATHS
unset_vars EXTERNAL
