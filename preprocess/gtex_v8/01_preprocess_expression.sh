#!/bin/bash -e

if [ "$PS1" ]; then echo -e "This script cannot be sourced. Use \"${BASH_SOURCE[0]}\" instead." ; return ; fi
CONFIGFILE=$1
if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ${BASH_SOURCE[0]} CONFIGFILE";
    exit 1;
fi
source ${CONFIGFILE}
source ${DATALOAD}
source ${EXTERNALLOAD}
source ../../main/PATHS

#---- Include functions
# source ${UTILSDIR}/gx_preproc_string
# source ${UTILSDIR}/submit_job

# Define all preprocessing scripts
PREPROC_SCRIPTDIR="${PWD}/scripts"
PREPROC_UTILSDIR="${PWD}/utils"
SELECTSAMPLEPY="${PREPROC_SCRIPTDIR}/select_samples_from_tissue.py"
GTEXNORMALIZEPY="${PREPROC_SCRIPTDIR}/gtex_normalization.py"  # do_expression_normalization.py
COMPILEAGECOVPY="${PREPROC_SCRIPTDIR}/compile_age_covariate.py"
CORRECTPY="${PREPROC_SCRIPTDIR}/correct_covariates.py"
PEERSCRIPT_R="${PREPROC_SCRIPTDIR}/PEER.R"
GENCODEFILTERPY="${PREPROC_SCRIPTDIR}/filter_gencode_expr.py"

PYPHASER="${PREPROC_SCRIPTDIR}/process_phASER_counts.py"
PYTPMS="${PREPROC_SCRIPTDIR}/calculated_TPMs.py"

# Get age covariate only once
mkdir -p $COVOUTDIR;
grep -v -P "^#" ${SRCSUBJCT} | cut -f 2,5,15 | grep -i gtex > ${AGE_COVARIATE_FILE}

# Define a random string for marking jobs of this batch
RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`

while IFS='' read -r LINE || [ -n "$LINE" ]; do
    if [[ $LINE =~ ^[^\#] ]]; then
        # Get the name of the tissue
        TFULL=$( echo "${LINE}" | cut -f 1 )
        TSHORT=$( echo "${LINE}" | cut -f 2 )
        TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )
        GXNPEERS=$( echo "${LINE}" | cut -f 3 )

        if [ "${USE_BSUB}" = "true" ]; then
            source utils/process_expression.sh
        else
            # Expression Files
            # NORMOUTFILE="${NORMOUTDIR}/${TSHORT}_normalized.txt"
            # LMOUTFILE="${LMOUTDIR}/${TSHORT}_lmcorrected.txt"
            # LMOUTFILE_AGE="${LMOUTDIR}/${TSHORT}_age_lmcorrected.txt"
            # Covariate Files
            COVARS="${COVOUTDIR}/${TSHORT}_nopeer_covariates.txt"
            COVARS_AGE="${COVOUTDIR}/${TSHORT}_nopeer_covariates_w_age.txt"
            COVARS_NOPC="${COVOUTDIR}/${TSHORT}_nopeer_covariates_nopc.txt"
            COVARS_AGE_NOPC="${COVOUTDIR}/${TSHORT}_nopeer_covariates_nopc_w_age.txt"

            if [ "${bSelectNormalize}" = "true" ]; then source ${PREPROC_UTILSDIR}/gx_preproc_01_select_normalize.sh; fi
            if [ "${bFormatCovariates}" = "true" ]; then source ${PREPROC_UTILSDIR}/gx_preproc_02_format_covariates.sh; fi
            if [ "${bLMcorrect}" = "true" ]; then source ${PREPROC_UTILSDIR}/gx_preproc_03_lmcorrect.sh; fi 
            if [ "${bPeerCorrect}" = "true" ]; then source ${PREPROC_UTILSDIR}/gx_preproc_04_peer.sh; fi 
            if [ "${bGencodeFilter}" = "true" ]; then source ${PREPROC_UTILSDIR}/gx_preproc_05_gencode_filter.sh; fi
        fi
    fi
done < ${TISSUEFILE}
