#!/bin/bash -e
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
# source ${UTILSDIR}/gx_preproc_string
# source ${UTILSDIR}/submit_job

# Define all preprocessing scripts
PREPROC_SCRIPTDIR="${PWD}/scripts"
PREPROC_UTILSDIR="${PWD}/utils"

PYTPMS="${PREPROC_SCRIPTDIR}/calculate_TPMs.py"
GTEXNORMALIZEPY="${PREPROC_SCRIPTDIR}/gtex_normalization.py"  # do_expression_normalization.py
GENCODEFILTERPY="${PREPROC_SCRIPTDIR}/filter_gencode_expr.py"

# Define a random string for marking jobs of this batch
RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`

SRCREAD="${GXDIR}/Counts.reproc.462.allgenes.txt"
OUTDIR="${GXOUTDIR}/tpms"
if [ ! -e $OUTDIR ]; then
    mkdir $OUTDIR
fi
OUTTPM="${GXOUTDIR}/tpms/TPMs_GEUVADIS.reproc.462.allgenes.txt"

echo ${PYENV} ${PYTPMS} --in ${SRCREAD} --out ${OUTTPM} --gtf ${GTFFILE}
# ${PYENV} ${PYTPMS} --in ${SRCREAD} --out ${OUTTPM} --gtf ${GTFFILE}


# Normalize
echo ${PYENV} ${GTEXNORMALIZEPY} --rpkm $OUTTPM --counts $SRCREAD --outdir $GXOUTDIR
# ${PYENV} ${GTEXNORMALIZEPY} --rpkm $OUTTPM --counts $SRCREAD --outdir $GXOUTDIR

# Filter protein coding and lncRNA
ALLGXFILES=()

EXPR_TYPES="tmm tpms"

# normalized uncorrected
ALLGXFILES+=( )

for EXPR_TYPE in ${EXPR_TYPES}; do
    if [ "${EXPR_TYPE}" == "tpms" ] || [ "${EXPR_TYPE}" == "rpkms" ]; then
        ALLGXFILES+=( "${GXOUTDIR}/${EXPR_TYPE}/${EXPR_TYPE}_qcfilter.txt" )
    else
        ALLGXFILES+=( "${GXOUTDIR}/${EXPR_TYPE}/${EXPR_TYPE}.txt" )
    fi
done

for THIS_GXFILE in ${ALLGXFILES[@]}; do
    echo $THIS_GXFILE
    if [ -e $THIS_GXFILE ]; then 
        echo ${PYENV} ${GENCODEFILTERPY} --gx ${THIS_GXFILE} \
                                    --dataset gtex \
                                    --gtf ${GENCODEFILE} \
                                    --biotype ${GXSELECTION} # --out ${TISSUEGXFILE} 
    fi
done