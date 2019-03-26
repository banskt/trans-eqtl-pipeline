#!/bin/bash -e

source DATA
source CONFIG
source EXTERNAL
source PATHS
source ${UTILSDIR}/submit_job

mkdir -p $RPKMOUTDIR;
mkdir -p $COVOUTDIR;
mkdir -p $PEEROUTDIR;

# Get age covariate only once
grep -v -P "^#" ${SRCSUBJCT} | cut -f 2,5,13 |grep -i gtex > $AGE_COVARIATE


while IFS='' read -r LINE || [ -n "$LINE" ]; do
    if [[ $LINE =~ ^[^\#] ]]; then
        # Get the name of the tissue
        TFULL=$( echo "${LINE}" | cut -f 1 )
        TSHORT=$( echo "${LINE}" | cut -f 2 )
        TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

        if [ "${USE_BSUB}" = "true" ]; then
            source ${UTILSDIR}/process_gxpression.sh
        else
            # Expression Files
            EXPRFILE="${NORMOUTDIR}/${TSHORT}_normalized.txt"
            LMOUTFILE="${LMOUTDIR}/${TSHORT}_lmcorrected.txt"
            LMOUTFILE_AGE="${LMOUTDIR}/${TSHORT}_age_lmcorrected.txt"

            # Covariate Files
            COVARS="${COVOUTDIR}/${TSHORT}_nopeer_covariates.txt"
            COVARS_AGE="${COVOUTDIR}/${TSHORT}_nopeer_covariates_w_age.txt"

            if [ "${bProcessRPKMs_and_normalize}" = "true" ]; then source ${UTILSDIR}/process_rpkms_and_normalize.sh; fi
            if [ "${bformatCovariates}" = "true" ]; then source ${UTILSDIR}/format_covariates.sh; fi
            if [ "${bLMcorrect}" = "true" ]; then source ${UTILSDIR}/lm_correct.sh; fi
            if [ "${GTEX_PEER_CORRECTION}" = "true" ]; then
                if [ "${CORRECT_AGE}" = "true" ]; then
                    source ${UTILSDIR}/peer_correction.sh $EXPRFILE $TSHORT "${PEEROUTDIR}_w_age" $COVARS_AGE;
                else
                    source ${UTILSDIR}/peer_correction.sh $EXPRFILE $TSHORT $PEEROUTDIR $COVARS;
                fi
            fi
            if [ "${GENCODE_FILTER}" = "true" ]; then source ${UTILSDIR}/gencode_filter.sh; fi
        fi
    fi
done < ${TISSUEFILE}
