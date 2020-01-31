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
PYDHS="${SCRIPTDIR}/check_dhs_enrichment.py"

for EXPR_CORR in ${EXPRESSIONS}; do
    for MDATA in ${DATASETS}; do
        for ANNOT_TYPE in ${ANNOTS}; do 
            source ${DATALOAD}

            JOBDIR_DATA="${JOBSUBDIR}/dhs_enrichment/${MDATA}"
            if [ -d ${JOBDIR_DATA} ]; then rm -rf ${JOBDIR_DATA}; fi; mkdir -p ${JOBDIR_DATA}
            JOBNAME="dhs_${DATATYPE}_${RANDSTRING}"

            OUTDIR_DATA="${OUTDIR}/${MDATA}"
            DESC="${EXPR_CORR}_${ANNOT_TYPE}"
            if [ $KNN -gt 0 ]; then
                DESC="${DESC}_knn${KNN}"
            fi

            echo "${DESC}"
            if [ ${DYNAMIC_SB} == "true" ]; then
                for KEFF in $KEFFS; do
                    if [ "${bTejaas}" = "true" ] && [ "${bTjsRandom}" = "true" ];   then
                        RESULTS="${OUTDIR_DATA}/tejaas/permnull_sb${SBETA}/chr{:d}/rr.txt.ld_prune"
                        NULL_RESULTS="${OUTDIR_DATA}/tejaas_rand/permnull_sb${SBETA}/chr{:d}/rr.txt.ld_prune" 
                    fi
                
                    ANNOTSNPS="${OUTDIR_DATA}/tejaas/permnull_sb${SBETA}/SNPs_annots.txt"
                    PLOT_OUTFILE="${OUTDIR}/../${MDATA}_tejaas_${SBETA}_${DESC}.png"

                    ${PYENV} ${PYDHS} --snpfile ${RESULTS} \
                                    --dhsfile ${DHS_FILE}\
                                    --annotfile ${ANNOTSNPS} \
                                    --outfile ${PLOT_OUTFILE} \
                                    --desc ${DESC} \
                                    --null-snpfile ${NULL_RESULTS}

                done
            else
                for SBETA in $TEJAAS_SIGMA_BETA_PERM; do
                    if [ "${bTejaas}" = "true" ] && [ "${bTjsRandom}" = "true" ];   then
                        RESULTS="${OUTDIR_DATA}/tejaas/permnull_sb${SBETA}/chr{:d}/rr.txt.ld_prune"
                        NULL_RESULTS="${OUTDIR_DATA}/tejaas_rand/permnull_sb${SBETA}/chr{:d}/rr.txt.ld_prune" 
                    fi
                
                    ANNOTSNPS="${OUTDIR_DATA}/tejaas/permnull_sb${SBETA}/SNPs_annots.txt"
                    PLOT_OUTFILE="${OUTDIR}/../${MDATA}_tejaas_${SBETA}_${DESC}.png"

                    ${PYENV} ${PYDHS} --snpfile ${RESULTS} \
                                    --dhsfile ${DHS_FILE}\
                                    --annotfile ${ANNOTSNPS} \
                                    --outfile ${PLOT_OUTFILE} \
                                    --desc ${DESC} \
                                    --null-snpfile ${NULL_RESULTS}

                done
            fi
        done
    done
done

