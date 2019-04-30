#!/bin/bash

## select which files to be filtered

GXSELECTION_STRING=$( echo ${GXSELECTION} | sed 's/ /_/g' )
GXFILENAME_FMT_THIS="${GXFILENAME_FMT/\[SELECTION\]/${GXSELECTION_STRING}}"
ALLGXFILES=()
PREPROC_STRINGS=()

if [ ! "${bPeerCorrect}" = "true" ]; then
    ## No covariate correction
    THISGXFILE="${NORMOUTFILE}"
    PREPROC_STRING=$( gx_preproc_string true false false 0 )
    ALLGXFILES+=( "$THISGXFILE" )
    PREPROC_STRINGS+=( "$PREPROC_STRING" )

    ## With covariate correction
    for bAgeCorrect in ${GXAGECORR}; do
        if [ "${bAgeCorrect}" = "true" ]; then
            THISGXFILE="${LMOUTFILE_AGE}"
        else 
            THISGXFILE="${LMOUTFILE}"
        fi
        for NPEER in ${GXNPEERS}; do
            PREPROC_STRING=$( gx_preproc_string true true ${bAgeCorrect} 0 )
            PREPROC_STRINGS+=( "$PREPROC_STRING" )
            ALLGXFILES+=( "$THISGXFILE" )
        done
    done

else
    ## No covariate correction
    for NPEER in ${GXNPEERS}; do
        THISGXFILE="${PEEROUTDIR}/normalized/${TSHORT}_${NPEER}_PEER_residuals.txt"
        PREPROC_STRING=$( gx_preproc_string true false false ${NPEER} )
        PREPROC_STRINGS+=( "$PREPROC_STRING" )
        ALLGXFILES+=( "$THISGXFILE" )
    done

    ## With covariate correction
    for bAgeCorrect in ${GXAGECORR}; do
        if [ "${bAgeCorrect}" = "true" ]; then
            THIS_PEEROUTDIR="${PEEROUTDIR}/covar_withage"
        else 
            THIS_PEEROUTDIR="${PEEROUTDIR}/covar"
        fi
        for NPEER in ${GXNPEERS}; do
            THISGXFILE="${THIS_PEEROUTDIR}/${TSHORT}_${NPEER}_PEER_residuals.txt"
            PREPROC_STRING=$( gx_preproc_string true true ${bAgeCorrect} ${NPEER} )
            PREPROC_STRINGS+=( "$PREPROC_STRING" )
            ALLGXFILES+=( "$THISGXFILE" )
        done
    done
fi

INDEX=0
for PREPROC in ${PREPROC_STRINGS[@]}; do
    THIS_PREPROC_FMT="${GXFILENAME_FMT_THIS/\[PREPROC\]/${PREPROC}}"
    TISSUEGXFILE="${THIS_PREPROC_FMT/\[TISSUE\]/${TSHORT}}"
    ${PYENV} ${GENCODEFILTERPY} --gx ${ALLGXFILES[${INDEX}]} --donors ${DONORFILE} --dataset gtex --out ${TISSUEGXFILE} --gtf ${GENCODEFILE} --biotype ${GXSELECTION}
    INDEX=$(( INDEX + 1 ))
done
