#!/bin/bash

source ${UTILSDIR}/submit_job

THISJOBDEPS="None" # no job dependencies
SPECIFIC_JOBSUBDIR="${PREPROCDIR}/framingham/jobsubs"
SPECIFIC_OUTDIR="${SRCDIR}/genotypes/merged_dosages"

INFILE="${SRCDIR}/genotypes/[CNST]/chr[CHRM]_[CNST].fhs.dosages.txt.gz"
INFILE_JOE="${SRCDIR}/genotypes/[CNST]/chr[CHRM]_[CNST].fhs.dosages.joehannes.txt.gz"
SAMPLES="${SRCDIR}/genotypes/[CNST]/chr[CHRM]_[CNST].fhs.dosages.sample"

# if [ -d ${SPECIFIC_JOBSUBDIR} ]; then rm -rf ${SPECIFIC_JOBSUBDIR}; fi; mkdir -p ${SPECIFIC_JOBSUBDIR}

if [ -d ${SPECIFIC_OUTDIR} ]; then
    echo "Directory exists!: ${SPECIFIC_OUTDIR}"
    exit
else
    mkdir -p $SPECIFIC_OUTDIR
fi

OUTFILE="${SPECIFIC_OUTDIR}/chr[CHRM].fhs.dosages.txt.gz"
OUTFILE_JOE="${SPECIFIC_OUTDIR}/chr[CHRM].fhs.dosages.joehannes.txt.gz"
OUTFILE_SAMPLES="${SPECIFIC_OUTDIR}/chr[CHRM].fhs.dosages.sample"

CONSENT_GROUPS="c1 c2"
echo "Overwriting options for BOTH consent groups C1 and C2"

for CHRM in ${CHRNUMS}; do
    COMMAND="paste -d \" \""
    COMMAND_JOE="paste -d \" \""
    CAT_CMMD="cat"

    GENOTYPE=${INFILE/\[CHRM\]/${CHRM}}
    GENOTYPE_JOE=${INFILE_JOE/\[CHRM\]/${CHRM}}
    CHRM_OUTFILE=${OUTFILE/\[CHRM\]/${CHRM}}
    CHRM_OUTFILE_JOE=${OUTFILE_JOE/\[CHRM\]/${CHRM}}

    SAMPLE_CHRM=${SAMPLES/\[CHRM\]/${CHRM}}
    OUTSAMPLE_CHRM=${OUTFILE_SAMPLES/\[CHRM\]/${CHRM}}

    JOBNAME="fhs_gt_merge_${CHRM}_${RANDSTRING}"
    for CONSENT in ${CONSENT_GROUPS}; do
        GT_CHRM_CNST=${GENOTYPE//\[CNST\]/${CONSENT}}
        GT_CHRM_CNST_JOE=${GENOTYPE_JOE//\[CNST\]/${CONSENT}}
        SAMPLE_CHRM_CNST=${SAMPLE_CHRM//\[CNST\]/${CONSENT}}

        if [ $CONSENT = "c1" ]; then 
            COMMAND="${COMMAND} <(gzip -cd $GT_CHRM_CNST)"
            COMMAND_JOE="${COMMAND_JOE} <(gzip -cd $GT_CHRM_CNST_JOE)"
            CAT_CMMD="${CAT_CMMD} $SAMPLE_CHRM_CNST"
        else  
            COMMAND="${COMMAND} <(gzip -cd $GT_CHRM_CNST | cut -d \" \" -f 7-)"
            COMMAND_JOE="${COMMAND_JOE} <(gzip -cd $GT_CHRM_CNST_JOE | cut -d \" \" -f 7-)"
            CAT_CMMD="${CAT_CMMD} <(tail -n +3 $SAMPLE_CHRM_CNST )"
        fi
    done

    CAT_CMMD="$CAT_CMMD > $OUTSAMPLE_CHRM"
    COMMAND="$COMMAND |gzip > $CHRM_OUTFILE"
    COMMAND_JOE="$COMMAND_JOE |gzip > $CHRM_OUTFILE_JOE"

    # echo $CAT_CMMD
    # echo $COMMAND
    # echo $COMMAND_JOE

    sed -e "s@_JOB_NAME@${JOBNAME}@g;
            s@_CAT_@${CAT_CMMD}@g;
            s@_CMD_@${COMMAND}@g;
            s@_CMMD_JOE_@${COMMAND_JOE}@g;
            " ${MASTER_BSUBDIR}/fhs_gt_merge.bsub > ${SPECIFIC_JOBSUBDIR}/${JOBNAME}.bsub

        submit_job ${SPECIFIC_JOBSUBDIR} ${JOBNAME} ${THISJOBDEPS}
done