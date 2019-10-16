#!/bin/bash

if [ "$PS1" ]; then echo -e "This script cannot be sourced. Use \"${BASH_SOURCE[0]}\" instead." ; return ; fi

CONFIGFILE=$1
if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ${BASH_SOURCE[0]} CONFIGFILE";
    exit 1;
fi

source ${CONFIGFILE}
source ${EXTERNALLOAD}
source PATHS
source ${UTILSDIR}/submit_job
source ${UTILSDIR}/unset_vars

RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
THISJOBDEPS="None"

for MDATA in ${LD_DATASETS}; do
    source ${DATALOAD}
    OUTDIR_DATA="${GENO_DIR}/ldmap_${LDWINDOW}_${LDMIN_R2}"
    if [ ! -d ${OUTDIR_DATA} ]; then mkdir -p ${OUTDIR_DATA}; fi

    JOBDIR_DATA="${JOBSUBDIR}/ldmap_${LDWINDOW}_${LDMIN_R2}/${DATATYPE}"
    if [ -d ${JOBDIR_DATA} ]; then rm -rf ${JOBDIR_DATA}; fi; mkdir -p ${JOBDIR_DATA}

    if [ -f ${LD_REGIONS_FILE} ]; then
        while IFS="" read -r line || [ -n "$p" ]
        do
            printf '%s\n' "$line"
            CHRM=`echo ${line} | cut -d" " -f 1`
            FROM=`echo ${line} | cut -d" " -f 2`
            TO=`echo ${line} | cut -d" " -f 3`

            JOBNAME="ld_${DATATYPE}_${CHRM}_${FROM}_${TO}_${RANDSTRING}"
            VCFFILE="${GENO_VCF_FMT/\[CHRM\]/${CHRM}}"
            LDOUTFILE="${OUTDIR_DATA}/chr${CHRM}_${FROM}_${TO}_${MDATA} "

            sed -e "s|_JOB_NAME|${JOBNAME}|g;
                    s|_VCFTOOL_|${VCFTOOLS}|g;
                    s|_GT_FILE_|${VCFFILE}|g;
                    s|_WINDOW_|${LDWINDOW}|g;
                    s|_MIN_R2_|${LDMIN_R2}|g;
                    s|_CHRM_|${CHRM}|g;
                    s|_OUTF_|${LDOUTFILE}|g;
                    s|_KEEP_|${LDKEEPFILE}|g;
                    s|_FROM_|${FROM}|g;
                    s|_TO_|${TO}|g;
                   " ${MASTER_BSUBDIR}/create_ldmap_regions.bsub > ${JOBDIR_DATA}/${JOBNAME}.bsub
            
            # bsub < ${JOBNAME}.bsub
            submit_job ${JOBDIR_DATA} ${JOBNAME} ${THISJOBDEPS}
        done < $LD_REGIONS_FILE
    else
        for CHRM in ${CHRNUMS}; do

            JOBNAME="ld_${DATATYPE}_${CHRM}_${RANDSTRING}"
            VCFFILE="${GENO_VCF_FMT/\[CHRM\]/${CHRM}}"
            LDOUTFILE="${OUTDIR_DATA}/chr${CHRM}_${MDATA} "

            sed -e "s|_JOB_NAME|${JOBNAME}|g;
                    s|_VCFTOOL_|${VCFTOOLS}|g;
                    s|_GT_FILE_|${VCFFILE}|g;
                    s|_WINDOW_|${LDWINDOW}|g;
                    s|_MIN_R2_|${LDMIN_R2}|g;
                    s|_OUTF_|${LDOUTFILE}|g;
                    s|_KEEP_|${LDKEEPFILE}|g;
                   " ${MASTER_BSUBDIR}/create_ldmap.bsub > ${JOBDIR_DATA}/${JOBNAME}.bsub
            
            # bsub < ${JOBNAME}.bsub
            submit_job ${JOBDIR_DATA} ${JOBNAME} ${THISJOBDEPS}
        done
    fi
    unset_vars ${DATALOAD}
done
