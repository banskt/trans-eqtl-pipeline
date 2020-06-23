#!/bin/bash -e

CHRMNUMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
# CHRMNUMS="9"

source "${PWD}/utils/submit_job"
VCFTOOLS="/usr/users/fsimone/bin/vcftools_bins/bin/vcftools"

GTEXDIR_V8="/cbscratch/franco/datasets/gtex_v8"
INFILE_FMT="${GTEXDIR_V8}/genotypes/vcfs_SHAPEIT2/0.01/GTEX_v8_2020-02-21_WGS_838Indiv_Freeze.SHAPEIT2_phased_NoMissingGT_SNPfilter_MAF0.01_chr[CHRM].vcf.gz"
OUTDIR="${GTEXDIR_V8}/genotypes/vcfs_SHAPEIT2/fst"
if [ ! -e $OUTDIR ]; then
    mkdir -p $OUTDIR
fi

THISJOBDEPS="None" # no job dependencies
SPECIFIC_JOBSUBDIR="${PWD}/jobsubs/fst"
if [ -d ${SPECIFIC_JOBSUBDIR} ]; then rm -rf ${SPECIFIC_JOBSUBDIR}; fi; mkdir -p ${SPECIFIC_JOBSUBDIR}

# Define a random string for marking jobs of this batch
RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`

for CHRM in $CHRMNUMS; do 

    POP1="${GTEXDIR_V8}/genotypes/gtex_v8_eur.list"
    POP2="${GTEXDIR_V8}/genotypes/gtex_v8_afr.list"
    POP3="${GTEXDIR_V8}/genotypes/gtex_v8_oth.list"

    CHRM_IN=${INFILE_FMT/\[CHRM\]/${CHRM}}

    #### 3 populations ####
    JOBNAME="gtex_v8_fst_3pop_[CHRM]_${RANDSTRING}"
    OUTFILE_FMT="GTEx_v8_SHAPEIT2_EUR-AFR-OTH-chr[CHRM]"
    CHRM_OUT=${OUTFILE_FMT/\[CHRM\]/${CHRM}}
    CHRM_JOBNAME=${JOBNAME/\[CHRM\]/${CHRM}}

    OUTFILE="${OUTDIR}/${CHRM_OUT}"

    sed -e "s|_VCFTOOL_|${VCFTOOLS}|g;
            s|_GTFILE_|${CHRM_IN}|g;
            s|_JOB_NAME|${CHRM_JOBNAME}|g;
            s|_OUTFILE_|${OUTFILE}|g;
            s|_POP1_|${POP1}|g;
            s|_POP2_|${POP2}|g;
            s|_POP3_|${POP3}|g
            " ${PWD}/bsubfiles/calc_fst.slurm > ${SPECIFIC_JOBSUBDIR}/${CHRM_JOBNAME}.slurm

    submit_job ${SPECIFIC_JOBSUBDIR} ${CHRM_JOBNAME} ${THISJOBDEPS}

    #### 2 populations ####
    JOBNAME="gtex_v8_fst_2pop_[CHRM]_${RANDSTRING}"
    OUTFILE_FMT="GTEx_v8_SHAPEIT2_EUR-AFR-chr[CHRM]"
    CHRM_OUT=${OUTFILE_FMT/\[CHRM\]/${CHRM}}
    CHRM_JOBNAME=${JOBNAME/\[CHRM\]/${CHRM}}

    OUTFILE="${OUTDIR}/${CHRM_OUT}"

    POP3=""
    sed -e "s|_VCFTOOL_|${VCFTOOLS}|g;
            s|_GTFILE_|${CHRM_IN}|g;
            s|_JOB_NAME|${CHRM_JOBNAME}|g;
            s|_OUTFILE_|${OUTFILE}|g;
            s|_POP1_|${POP1}|g;
            s|_POP2_|${POP2}|g;
            s|_POP3_|${POP3}|g
            " ${PWD}/bsubfiles/calc_fst.slurm > ${SPECIFIC_JOBSUBDIR}/${CHRM_JOBNAME}.slurm

    submit_job ${SPECIFIC_JOBSUBDIR} ${CHRM_JOBNAME} ${THISJOBDEPS}

done;


# for chrm in {1..22}; do tail -n +2 GTEx_v8_SHAPEIT2_EUR-AFR-chr${chrm}.weir.fst >> GTEx_EUR_AFR.weir.fst; done
# cat GTEx_EUR_AFR.weir.fst | sed 's/chr//g' > GTEx_EUR_AFR_fix.weir.fst 