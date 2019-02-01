#!/bin/bash

CONFIGFILE=$1

#---- Run this code with CONFIGFILE -----#
if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ./02_preprocess_genotype.sh CONFIGFILE"
    exit 1
fi

#---- It also needs the DATA file, which is hard-coded as of now.
#---- However, the original data is always the same
#---- and therefore, unlike CONFIG, several options are not needed.
#---- Finally it needs the path or directory structure of the pipeline
source ${CONFIGFILE}
source DATA
source EXTERNAL
source ../../main/PATHS

#---- Include functions
source ${UTILSDIR}/gx_preproc_string
source ${UTILSDIR}/submit_job
source ${UTILSDIR}/add_deps
source utils/process_gx_rows

VCF_SELECTSORT_PY="${SCRIPTDIR}/vcf_select_sort_samples.py"

PREPROC_STRING=$( gx_preproc_string ${GXNORMALIZE[0]} ${GXLMCORR[0]} ${GXNPEER[0]} )
GXFILENAME_FMT_THIS="${GXFILENAME_FMT/\[SELECTION\]/${GXSELECTION[0]}}"
GXFILENAME_FMT_THIS="${GXFILENAME_FMT_THIS/\[PREPROC\]/${PREPROC_STRING}}" # gene expression filename

JOBDEPS="None" # used for controlling job dependencies
RANDSTRING=$( cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1 )
INFILE_BASE="${GTOUTDIR}/all_samples/${GTFILE_BASENAME}_maf${MAFMIN#*.}"
if [ "${REMOVE_INDELS}" = "true" ]; then INFILE_BASE="${INFILE_BASE}_noindels"; fi
if [ "${REMOVE_AMBIGUOUS}" = "true" ]; then INFILE_BASE="${INFILE_BASE}_noambig"; fi


#---- Create genotypes for each tissue separately (!), 
#---- specifically including only those samples 
#---- which are present in the gene expression of those tissues
#---- and in the same order as in the expression file.
#---- Dosage format is required for TEJAAS and MatrixEQTL
#---- Bed format with no missing genotype is required for GNetLMM.

while IFS='' read -r LINE || [ -n "${LINE}" ]; do

    # Get the name of the tissue
    TFULL=$( echo "${LINE}" | cut -f 1 )
    TSHORT=$( echo "${LINE}" | cut -f 2 )
    TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

    # Which expression file?
    EXPRFILE="${GXFILENAME_FMT_THIS/\[TISSUE\]/${TSHORT}}"
    GTOUT_THIS="${GTOUTDIR}/${TBASE}"

    # Proceed only if the expression file exists
    if [ -e ${EXPRFILE} ] ; then

        # keep only common samples, sorted according to expression file
        # (.bed, .bim and .fam files)
        # gene expression is converted to GNetLMM format
        # (.matrix, .cols and .rows files)
        echo "### TISSUE: ${TFULL} ####"
        source utils/preprocess_gnetlmm.sh
    fi
    
done < ${TISSUEFILE}
