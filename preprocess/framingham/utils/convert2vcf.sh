#!/bin/bash
source ${UTILSDIR}/submit_job

THISJOBDEPS="None" # no job dependencies
SPECIFIC_JOBSUBDIR="${PREPROCDIR}/framingham/jobsubs"
SPECIFIC_INDIR="${SRCDIR}/genotypes/merged_dosages"
SPECIFIC_OUTDIR="${SRCDIR}/genotypes/vcfs"
if [ -d ${SPECIFIC_JOBSUBDIR} ]; then rm -rf ${SPECIFIC_JOBSUBDIR}; fi; mkdir -p ${SPECIFIC_JOBSUBDIR}

# if [ -d ${SPECIFIC_OUTDIR} ]; then
#     echo "Directory exists!: ${SPECIFIC_OUTDIR}"
#     exit
# else
#     mkdir -p $SPECIFIC_OUTDIR
# fi

INFILE="${SPECIFIC_INDIR}/chr[CHRM].fhs.dosages.txt.gz"
INFILE_JOE="${SPECIFIC_INDIR}/chr[CHRM].fhs.dosages.joehannes.txt.gz"
INFILE_SAMPLES="${SPECIFIC_INDIR}/chr[CHRM].fhs.dosages.sample"

OUTVCF="${SPECIFIC_OUTDIR}/chr[CHRM].fhs.vcf.gz"
OUTVCF_JOE="${SPECIFIC_OUTDIR}/chr[CHRM].fhs.joehannes.vcf.gz"

for CHRM in ${CHRNUMS}; do
    JOBNAME="fhs_gt_vcf_${CHRM}_${RANDSTRING}"

    CHRM_INFILE=${INFILE/\[CHRM\]/${CHRM}}
    CHRM_INFILE_JOE=${INFILE_JOE/\[CHRM\]/${CHRM}}

    SAMPLE_CHRM=${INFILE_SAMPLES/\[CHRM\]/${CHRM}}

    CHRM_OUT=${OUTVCF/\[CHRM\]/${CHRM}}
    CHRM_OUT_JOE=${OUTVCF_JOE/\[CHRM\]/${CHRM}}

    sed -e "s|_JOB_NAME|${JOBNAME}|g;
            s|_PYENV_|${PYENV}|g;
            s|_SMPLE_|${SAMPLE_CHRM}|g;
            s|_CALLT_|${VCF_CALLTHRES}|g;
            s|_CONVERT_|${CONVERTPY}|g;
            s|_NOTDOSE_|${NOTDOSAGE}|g;
            s|_INFILE_|${CHRM_INFILE}|g;
            s|_JOEFILE_|${CHRM_INFILE_JOE}|g;
            s|_OUTFILE_|${CHRM_OUT}|g;
            s|_OUTJOE_|${CHRM_OUT_JOE}|g;
            " ${MASTER_BSUBDIR}/fhs_convert2vcf.bsub > ${SPECIFIC_JOBSUBDIR}/${JOBNAME}.bsub

    submit_job ${SPECIFIC_JOBSUBDIR} ${JOBNAME} ${THISJOBDEPS}

done;

