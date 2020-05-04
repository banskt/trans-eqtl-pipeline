#!/bin/bash
#SBATCH -p medium
#SBATCH --qos=short
#SBATCH -t 0-2:00:00
#SBATCH -n 1
#SBATCH -C scratch
#SBATCH -J _JOB_NAME
#SBATCH -o _JOB_NAME.out
#SBATCH -e _JOB_NAME.err
#SBATCH -A all

SPECIFIC_OUTDIR=_RES_DIR_
RESOUTDIR=_COL_DIR_
CHRM_NTOT_FILE=_CHR_NTT_
CHUNK_REDUCE_UTIL=_TJC_RED_
CHRNUMS=_CHR_NUM_
PP_TEJAAS_PVALCUT=_TJ_PCUT_
NMAX=_NMAX_JB_
bChunkReduce=_bCHNK_R_
bPostResult=_bPST_RS_

RANDSTRING=$( cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1 )
source ${CHUNK_REDUCE_UTIL}

#================================
# Collect the results from chunks
#================================

if [ "${bChunkReduce}" = "true" ]; then
    for CHRM in ${CHRNUMS}; do
        CHRINDX=$((CHRM - 1))

        # Find number of SNPs in this chromosome
        NTOT=$( sed "${CHRM}q;d" ${CHRM_NTOT_FILE} | awk '{print $2}' )

        # Find number of jobs submitted
        NCHUNK=$( echo $(( (NTOT + NMAX - 1)/NMAX )) )
        echo "chr${CHRM} --> ${NCHUNK} chunks"

        tejaas_chunk_reduce "${SPECIFIC_OUTDIR}/chr${CHRM}" ${NCHUNK}
    done
fi

#================================
# Prepare trans-eQTL results
#================================

if [ "${bPostResult}" = "true" ]; then

    if [ ! -d ${RESOUTDIR} ]; then mkdir -p ${RESOUTDIR}; fi
    
    __RROUT="${SPECIFIC_OUTDIR}/trans_eqtls.txt"
    __TGOUT="${SPECIFIC_OUTDIR}/target_genes.txt"
    __TGKNN="${SPECIFIC_OUTDIR}/target_genes_knn.txt"

    echo -e "ID\tChrom\tPos\tMAF\tQ\tMu_Q\tSigma_Q\tP-Value" > ${__RROUT}
    echo -e "SNP_ID\tGene_ID\tP-Value" > ${__TGOUT}
    echo -e "SNP_ID\tGene_ID\tP-Value" > ${__TGKNN}
    echo -ne "Creating list of trans-eQTLs and target genes from "

    for CHRM in ${CHRNUMS}; do
        echo -ne "${CHRM} "
        __CHROUTDIR="${SPECIFIC_OUTDIR}/chr${CHRM}"
        __TMPDIR="${__CHROUTDIR}/tmp_${RANDSTRING}"

        # remove temporary directory if any
        if [ -d ${__TMPDIR} ]; then rm -rf ${__TMPDIR}; fi; mkdir -p ${__TMPDIR};

        # select trans-eQTLs below cutoff
        cat ${__CHROUTDIR}/rr.txt | awk -v pcut=${PP_TEJAAS_PVALCUT} '$8 <= pcut {print}' > ${__TMPDIR}/rr_chr${CHRM}.txt
        cat ${__TMPDIR}/rr_chr${CHRM}.txt >> ${__RROUT}

        # create a list of trans-eQTL IDs
        cut -f1 ${__TMPDIR}/rr_chr${CHRM}.txt > ${__TMPDIR}/teqtl_chr${CHRM}.txt

        # select the trans-eGenes for the above list of trans-eQTLs
        grep -Fwf ${__TMPDIR}/teqtl_chr${CHRM}.txt ${__CHROUTDIR}/gene_snp_list.txt     | awk '{printf "%s\t%s\t%g\n", $2, $1, $3}' >> ${__TGOUT}
        grep -Fwf ${__TMPDIR}/teqtl_chr${CHRM}.txt ${__CHROUTDIR}/gene_snp_list_knn.txt | awk '{printf "%s\t%s\t%g\n", $2, $1, $3}' >> ${__TGKNN}

        # remove the temporary direcotory
        rm -rf ${__TMPDIR}
    done
    echo "Done."

    cp ${__RROUT} ${RESOUTDIR}/
    cp ${__TGOUT} ${RESOUTDIR}/
    cp ${__TGKNN} ${RESOUTDIR}/

fi
