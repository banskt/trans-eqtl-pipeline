#!/bin/bash

THISJOBDEPS="None"
SPECIFIC_JOBSUBDIR="${JOBSUBDIR}/preprocess/gtex/${TBASE}"
if [ -d ${SPECIFIC_JOBSUBDIR} ]; then rm -rf ${SPECIFIC_JOBSUBDIR}; fi; mkdir -p ${SPECIFIC_JOBSUBDIR}

if [ "${CREATE_EXPRESSION}" = "true" ]; then

    GXOUTFILE="$( dirname ${EXPRFILE} )/$( basename ${EXPRFILE} .txt )"
    echo "Creating expression matrix ..."
    cut -f2- ${EXPRFILE} | tail -n +2 | sed 's/\t/\ /g' > ${GXOUTFILE}.matrix
    echo "Saving column information ..."
    echo "sample_id" > ${GXOUTFILE}.cols
    head -n 1 ${EXPRFILE} | cut -f 2- | tr "\t" "\n" >> ${GXOUTFILE}.cols
    echo "Saving row information and acquiring metadata from GENCODE."
    echo "This may take some time. Please wait ..."
    process_gx_rows ${EXPRFILE} ${GXOUTFILE} ${GENEPOSFILE}
fi


if [ "${CREATE_GENOTYPE}" = "true" ]; then

    echo "Creating list of donors for ${TFULL}"
    if [ ! -d ${GTOUT_THIS} ]; then mkdir -p ${GTOUT_THIS}; fi
    DONORFILE=${GTOUT_THIS}/${TSHORT}.samples
    head -n 1 ${EXPRFILE} | cut -f 2- | tr "\t" "\n" > ${DONORFILE}
    wc -l $DONORFILE | awk '{for (i = 1; i <= $1; i++) print "1.0"}' > ${GTOUT_THIS}/ones.txt

    for CHRM in {1..22}; do
        
        JOBNAME="gtex_vcf_select_sort_samples_${TSHORT}_chr${CHRM}_${RANDSTRING}"
        INFILE="${INFILE_BASE}_chr${CHRM}_WARNimputedmissing.vcf.gz"
        OUTFILE="${GTOUT_THIS}/$( echo "$( basename ${INFILE} .vcf.gz )"_${TSHORT} )"
    
        sed -e "s|_JOB_NAME|${JOBNAME}|g;
                s|_BGZIP___|${BGZIP}|g;
                s|_TABIX___|${TABIX}|g;
                s|_PLINK2__|${PLINK2}|g;
                s|_IN_FILE_|${INFILE}|g;
                s|_OT_FILE_|${OUTFILE}|g;
                s|_DONR_FL_|${DONORFILE}|g;
                s|_VCSS_PY_|${VCF_SELECTSORT_PY}|g;
               " ${MASTER_BSUBDIR}/preprocess_gnetlmm_genotype.bsub > ${SPECIFIC_JOBSUBDIR}/${JOBNAME}.bsub
    
        echo "Job file: ${SPECIFIC_JOBSUBDIR}/${JOBNAME}.bsub"
        submit_job ${SPECIFIC_JOBSUBDIR} ${JOBNAME} ${THISJOBDEPS}
        JOBDEPS=$( add_deps "${JOBDEPS}" ${JOBNAME} )
    done
fi
