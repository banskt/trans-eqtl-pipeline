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
source ${UTILSDIR}/gx_preproc_string
source ${UTILSDIR}/submit_job

# Define all preprocessing scripts
PREPROC_SCRIPTDIR="${PWD}/scripts"
PREPROC_UTILSDIR="${PWD}/utils"
SELECTSAMPLEPY="${PREPROC_SCRIPTDIR}/select_samples_from_tissue.py"
PREPROCPY="${PREPROC_SCRIPTDIR}/preprocess_expression.py"
COMPILEAGECOVPY="${PREPROC_SCRIPTDIR}/compile_age_covariate.py"
PEERSCRIPT_R="${PREPROC_SCRIPTDIR}/PEER.R"
GENCODEFILTERPY="${PREPROC_SCRIPTDIR}/filter_gencode_expr.py"

# Get age covariate only once
mkdir -p $COVOUTDIR;
grep -v -P "^#" ${SRCSUBJCT} | cut -f 2,5,13 | grep -i gtex > ${AGE_COVARIATE_FILE}

# Define a random string for marking jobs of this batch
RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`

while IFS='' read -r LINE || [ -n "$LINE" ]; do
    if [[ $LINE =~ ^[^\#] ]]; then
        # Get the name of the tissue
        TFULL=$( echo "${LINE}" | cut -f 1 )
        TSHORT=$( echo "${LINE}" | cut -f 2 )
        TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

        if [ "${USE_BSUB}" = "true" ]; then
            source utils/process_expression.sh
        else
            ## Create covariate files
            COVARS="${COVOUTDIR}/${TSHORT}_nopeer_covariates.txt"
            COVARS_AGE="${COVOUTDIR}/${TSHORT}_nopeer_covariates_w_age.txt"
            if [ "${bFormatCovariates}" = "true" ]; then
                echo "Processing Covariates for Tissue: $TFULL"
                # Selects only genotype PCs, sex and platform
                grep -v -i "inferred" ${GTEXCOV} > ${COVARS}
                ${PYENV} ${COMPILEAGECOVPY} --input ${COVARS} --age ${AGE_COVARIATE_FILE} --output ${COVARS_AGE}
            fi

            ## Select the tissue
            if [ "${bSelectTissue}" = "true" ]; then
                ${PYENV} ${SELECTSAMPLEPY} --input ${SRCTPM} --output ${TPMFILE} --tissue="$TFULL" --pheno ${SRCPHENO}
            fi

            ## Perform suggested normalizations
            if [ "${bNormalizeQC}" = "true" ]; then
                ${PYENV} ${PREPROCPY} -g ${TPMFILE} -c ${COVARS} -o ${GXPREFIX} -m ${PREPROC_METHODS}
            fi

            ## Perform PEER correction
            if [ "${bPeerCorrect}" = "true" ]; then
                for NPEER in ${NPEERCORR}; do
                    if [ ! "${NPEER}" = "0" ]; then
                        for PRMETHOD in ${PEERTRG_METHODS}; do
                            GXOUTDIR=$( dirname ${GXPRCFILE} )
                            GXBASENAME=$( basename ${GXPRCFILE} )
                            GXOUTPREFIX="${GXBASENAME%%.*}"
                
                            PEER_INFILE="${GXOUTDIR}/${GXOUTPREFIX}_${PRMETHOD}.txt"
                            PEER_OUTFILE="${GXOUTDIR}/${GXOUTPREFIX}_${PRMETHOD}_${NPEER}_PEER_residuals.txt"
                
                            Rscript ${PEERCRXN_R} ${PEER_INFILE} ${GXOUTPREFIX}_${PRMETHOD}_${NPEER} --n ${NPEER} -o ${GXOUTDIR}
                            ${MPYTHON} ${PREPROCPY} -g ${PEER_OUTFILE} -c ${CFOUTFILE} -o ${PEER_OUTFILE} -m 'norm'
                        done
                    fi
                done
            fi

            ## Filter using GENCODE
            if [ "${bGencodeFilter}" = "true" ]; then
                INDEX=0
                for PREPROC in ${PREPROC_STRINGS[@]}; do
                    THIS_PREPROC_FMT="${GXFILENAME_FMT_THIS/\[PREPROC\]/${PREPROC}}"
                    TISSUEGXFILE="${THIS_PREPROC_FMT/\[TISSUE\]/${TSHORT}}"
                    ${PYENV} ${GENCODEFILTERPY} --gx ${ALLGXFILES[${INDEX}]} --donors ${DONORFILE} --dataset gtex --out ${TISSUEGXFILE} --gtf ${GENCODEFILE} --biotype ${GXSELECTION}
                    INDEX=$(( INDEX + 1 ))
                done
            fi
        fi
    fi
done < ${TISSUEFILE}
