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

GTFILE_BASENAME=$( basename ${SRCVCF} ".vcf.gz" ) # gtex genotype file name
GTOUT_ALL="${GTOUTDIR}/all_samples" # output directory for all samples
if [ ! -d ${GTOUT_ALL} ]; then mkdir -p ${GTOUT_ALL}; fi

PREPROC_STRING=$( gx_preproc_string ${GXNORMALIZE[0]} ${GXLMCORR[0]} ${GXNPEER[0]} )
GXFILENAME_FMT_THIS="${GXFILENAME_FMT/\[SELECTION\]/${GXSELECTION[0]}}"
GXFILENAME_FMT_THIS="${GXFILENAME_FMT_THIS/\[PREPROC\]/${PREPROC_STRING}}" # gene expression filename

RANDSTRING=$( cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1 )
JOBDEPS="None" # used for controlling job dependencies
CONVERT_ANNOT="${SCRIPTDIR}/vcf_change_annot.py"
IMPUTE_MISSING="${SCRIPTDIR}/vcf_impute_missingGT.py"

#---- Run jobs for splitting genotype into chromosomes, filtering, changing header of VCF files and modifying annotations
#---- ${JOBDEPS} contain the job dependencies
if [ "${bSPlitVcf}" = "true" ]; then
    source utils/annotation_split.sh
    source utils/vcf_split_filter_headerchange_annotation.sh
fi

#---- Create genotypes for each tissue separately (!), 
#---- specifically including only those samples 
#---- which are present in the gene expression of those tissues
#---- and in the same order as in the expression file.
#---- Dosage format is required for TEJAAS and MatrixEQTL
#---- Bed format with no missing genotype is required for GNetLMM.

while IFS='' read -r LINE || [ -n "$line" ]; do

    # Get the name of the tissue
    TFULL=$( echo "${LINE}" | cut -f 1 )
    TSHORT=$( echo "${LINE}" | cut -f 2 )
    TBASE=$( echo ${TFULL} | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

    # Which expression file?
    EXPRFILE="${GXFILENAME_FMT_THIS/\[TISSUE\]/${TSHORT}}"
    GTOUT_THIS="${GTOUTDIR}/${TBASE}"
    if [ -e ${EXPRFILE} ] ; then # Proceed only if the expression file exists
        echo "Creating list of donors for ${TFULL}"
        if [ ! -d ${GTOUT_THIS} ]; then mkdir -p ${GTOUT_THIS}; fi
        DONORFILE=${GTOUT_THIS}/${TSHORT}.samples
        head -n 1 ${EXPRFILE} | cut -f 2- | tr "\t" "\n" > ${DONORFILE}
 
        # only common samples, sorted according to expression file
        # missing genotype imputed in PED for GNetLMM
        source utils/tissue_specific_genotype.sh

####        for CHRM in {21..22}; do
####            echo "${TFULL}: chromosome ${CHRM}"
####            INFILE="${GTOUT_ALL}/${GTFILE_BASENAME}_chr${CHRM}.vcf.gz"
####            VCFGZ_OUTFILE="${GTOUT_THIS}/${GTFILE_BASENAME}_chr${CHRM}_${TSHORT}.vcf.gz"
####            PLINK_OUTFILE="${GTOUT_THIS}/${GTFILE_BASENAME}_chr${CHRM}_${TSHORT}_nomissing"
####            ANNOTFILE="${GTOUT_ALL}/${GTFILE_BASENAME}_chr${CHRM}.annot"
####            # source utils/tissue_specific_genotype.sh
####            #echo "python select_sort_samples_vcf.py --input ${INFILE} --out ${VCFGZ_OUTFILE} --incl-samples ${DONORFILE} --annot ${ANNOTFILE}"
####            #echo "plink2 --vcf ${VCFGZ_OUTFILE} dosage=DS --make-bed --maf 0.1 --max-maf 0.9 --geno dosage --hard-call-threshold 0.499999 --out ${PLINK_OUTFILE} > ${PLINK_OUTFILE}_plink.log"
####        done
    fi
    
done < ${TISSUEFILE}
