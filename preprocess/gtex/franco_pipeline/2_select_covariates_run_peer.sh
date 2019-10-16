#!/bin/bash

source "./CONFIG.sh"

mkdir -p $COVDIR

while IFS='' read -r line || [ -n "$line" ]; do
    fullname=$(echo "$line" | cut -f 1 )
    shortname=$(echo "$line" | cut -f 2 )
    base=$(echo $fullname | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )
    GTEXCOV="${GTEXCOVDIR}/${base}_Analysis.v6p.covariates.txt"
    COVARS="${GTEXCOVDIR}/${shortname}_nopeer_covariates.txt"

    if [ -e $GTEXCOV ]; then
        echo "Processing Covariates for Tissue: $fullname"
        # Selects only genotype PCs, sex and platform
        grep -v -i "inferred" $GTEXCOV > $COVARS

        if [ $RUNPEER = true ]; then
            # for NCOV in `seq 5 5 $MAXNCOV`;
            for NCOV in `seq 0 1 1`;
            do
                # get PEER covariates correcting for covariates above
                PEERPREFIX="${shortname}_gtex.${NCOV}ncov"
                EXPRFILE="${EXPROUTDIR}/gtex.normalized.expression.${shortname}.txt"
                
                Rscript PEER.R $EXPRFILE $PEERPREFIX --n $NCOV --covar $COVARS -o $COVDIR
            done
        fi
    else
        echo "No GTEx covariate file for $fullname"
    fi
done < $TISSUEFILE