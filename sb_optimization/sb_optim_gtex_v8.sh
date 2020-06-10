#!/bin/bash

TISSUEFILE="/usr/users/fsimone/trans-eqtl-pipeline/main/tissues.txt"
OUTDIR="gtex_v8"


mkdir -p $OUTDIR
mkdir -p ${OUTDIR}_jobsubs

TSHORTS=""
declare -A TISSUEPEER
while IFS='' read -r LINE || [ -n "${LINE}" ]; do
    if [[ $LINE =~ ^[^\#] ]]; then
        # Get the name of the tissue
        TFULL=$( echo "${LINE}" | cut -f 1 )
        TSHORT=$( echo "${LINE}" | cut -f 2 )
        TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )
        GXNPEERS=$( echo "${LINE}" | cut -f 3 )
        TSHORTS="${TSHORTS} ${TSHORT}"
        TISSUEPEER[$TSHORT]="${GXNPEERS}"
    fi
done < ${TISSUEFILE}

for TISSUE in $TSHORTS; do
    sed "s|__TISSUE__|${TISSUE}|g;
         s|__OUTDIR__|${OUTDIR}|g;
        " script_gtex.slurm > "${OUTDIR}_jobsubs/sboptim_${TISSUE}.slurm"
    sbatch "${OUTDIR}_jobsubs/sboptim_${TISSUE}.slurm"
done
