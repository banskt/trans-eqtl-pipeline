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
source ${UTILSDIR}/submit_job
source ${UTILSDIR}/add_deps

GTOUT_ALL="${GTOUTDIR}/all_samples" # output directory for all samples
if [ ! -d ${GTOUT_ALL} ]; then mkdir -p ${GTOUT_ALL}; fi

JOBDEPS="None" # used for controlling job dependencies
CONVERT_ANNOT_PY="${SCRIPTDIR}/vcf_change_annot.py"
VCF_FILTER_PY="${SCRIPTDIR}/vcf_filter.py"
IMPUTE_MISSING_PY="${SCRIPTDIR}/vcf_impute_missingGT.py"

RANDSTRING=$( cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1 )
OUTFILE_BASE="${GTOUT_ALL}/${GTFILE_BASENAME}_maf${MAFMIN#*.}"
if [ "${REMOVE_INDELS}" = "true" ]; then OUTFILE_BASE="${OUTFILE_BASE}_noindels"; fi
if [ "${REMOVE_AMBIGUOUS}" = "true" ]; then OUTFILE_BASE="${OUTFILE_BASE}_noambig"; fi

#---- Run jobs for splitting genotype into chromosomes, filtering, changing header of VCF files and modifying annotations
#---- ${JOBDEPS} contain the job dependencies
source utils/annotation_split.sh
source utils/vcf_split_filter_headerchange_annotation.sh
