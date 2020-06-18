#!/bin/bash -e

CHRMNUMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
TISSUEFILE="/usr/users/fsimone/trans-eqtl-pipeline/main/tissues.txt"

TSHORTS=""
while IFS='' read -r LINE || [ -n "${LINE}" ]; do
    if [[ $LINE =~ ^[^\#] ]]; then
        # Get the name of the tissue
        TFULL=$( echo "${LINE}" | cut -f 1 )
        TSHORT=$( echo "${LINE}" | cut -f 2 )
        TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )
        GXNPEERS=$( echo "${LINE}" | cut -f 3 )
        TSHORTS="${TSHORTS} ${TSHORT}"
    fi
done < ${TISSUEFILE}


source "${PWD}/utils/submit_job"
VCFTOOLS="/usr/users/fsimone/bin/vcftools_bins/bin/vcftools"

GTEXDIR_V8="/cbscratch/franco/datasets/gtex_v8"
POP1="${GTEXDIR_V8}/genotypes/sample_files/gtex_v8_[TISSUE]_eur.list"
POP2="${GTEXDIR_V8}/genotypes/sample_files/gtex_v8_[TISSUE]_afr.list"
POP3="${GTEXDIR_V8}/genotypes/sample_files/gtex_v8_[TISSUE]_oth.list"
INFILE_FMT="${GTEXDIR_V8}/genotypes/vcfs_SHAPEIT2/0.01/GTEX_v8_2020-02-21_WGS_838Indiv_Freeze.SHAPEIT2_phased_NoMissingGT_SNPfilter_MAF0.01_chr[CHRM].vcf.gz"
OUTFILE_FMT="GTEx_v8_SHAPEIT2_EUR-AFR-chr[CHRM]"

THISJOBDEPS="None" # no job dependencies
SPECIFIC_JOBSUBDIR="${PWD}/jobsubs/fst"
if [ -d ${SPECIFIC_JOBSUBDIR} ]; then rm -rf ${SPECIFIC_JOBSUBDIR}; fi; mkdir -p ${SPECIFIC_JOBSUBDIR}

# Define a random string for marking jobs of this batch
RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
JOBNAME="gtex_v8_fst_[CHRM]_[TISSUE]_${RANDSTRING}"

for CHRM in $CHRMNUMS; do 
    for TISSUE in $TSHORTS; do
        CHRM_IN=${INFILE_FMT/\[CHRM\]/${CHRM}}
        CHRM_OUT=${OUTFILE_FMT/\[CHRM\]/${CHRM}}
        CHRM_JOBNAME=${JOBNAME/\[CHRM\]/${CHRM}}

        CHRM_JOBNAME=${CHRM_JOBNAME/\[TISSUE\]/${TISSUE}}

        POP1_T=${POP1/\[TISSUE\]/${TISSUE}}
        POP2_T=${POP2/\[TISSUE\]/${TISSUE}}
        POP3_T=${POP3/\[TISSUE\]/${TISSUE}}

        OUTDIR="${GTEXDIR_V8}/genotypes/vcfs_SHAPEIT2/fst/${TISSUE}"
        if [ ! -e $OUTDIR ]; then
            mkdir -p $OUTDIR
        fi

        OUTFILE="${OUTDIR}/${CHRM_OUT}"

        sed -e "s|_VCFTOOL_|${VCFTOOLS}|g;
                s|_GTFILE_|${CHRM_IN}|g;
                s|_JOB_NAME|${CHRM_JOBNAME}|g;
                s|_OUTFILE_|${OUTFILE}|g;
                s|_POP1_|${POP1_T}|g;
                s|_POP2_|${POP2_T}|g;
                s|_POP3_|${POP3_T}|g
                " ${PWD}/bsubfiles/calc_fst.slurm > ${SPECIFIC_JOBSUBDIR}/${CHRM_JOBNAME}.slurm

        submit_job ${SPECIFIC_JOBSUBDIR} ${CHRM_JOBNAME} ${THISJOBDEPS}
done;


# for chrm in {1..22}; do tail -n +2 GTEx_v8_SHAPEIT2_EUR-AFR-chr${chrm}.weir.fst >> GTEx_EUR_AFR.weir.fst; done
# cat GTEx_EUR_AFR.weir.fst | sed 's/chr//g' > GTEx_EUR_AFR_fix.weir.fst 