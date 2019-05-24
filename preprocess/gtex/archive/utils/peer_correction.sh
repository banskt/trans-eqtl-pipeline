#!/bin/bash -e

# This file can only be sourced

EXPRFILE=$1
PEERPREFIX=$2
THISOUTDIR=$3
COVARS=$4

mkdir -p $THISOUTDIR;

for NCOV in $NCOVS;
do
    # get PEER covariates correcting for covariates above   
    Rscript ${SCRIPTDIR}/PEER.R $EXPRFILE "${PEERPREFIX}_${NCOV}" --n $NCOV --covar $COVARS -o $THISOUTDIR
done

