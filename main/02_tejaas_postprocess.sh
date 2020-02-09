#!/bin/bash

while IFS='' read -r LINE || [ -n "$LINE" ]; do
    if [[ $LINE =~ ^[^\#] ]]; then

        # Get the name of the tissue
        TFULL=$( echo "${LINE}" | cut -f 1 )
        TSHORT=$( echo "${LINE}" | cut -f 2 )
        TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

        echo "${TSHORT} (${TBASE})"

        JOBNAME="postprocess_${TSHORT}"

        sed -e "s|_JOB_NAME|${JOBNAME}|g;
                s|__TBASE__|\"${TBASE}\"|g;
                s|__TSHORT__|\"${TSHORT}\"|g;
               " process_chunks_master.sbatch > ${JOBNAME}.sbatch

        #bsub < ${JOBNAME}.bsub
        #sbatch ${JOBNAME}.sbatch
        #sleep 3

    fi
done < tissue_table.txt
