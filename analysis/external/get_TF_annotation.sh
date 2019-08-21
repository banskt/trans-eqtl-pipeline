#!/bin/bash
TFCSV=$1
GTFFILE=$2
OUTFILE=$3

RANDSTRING=$( cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1 )
TMP_TF_FILE="tmp_tf_${RANDSTRING}.txt"
TMP_GTF_UNGZ="tmp_gencode_${RANDSTRING}.gtf"
TMP_SEL_FILE="tmp_sel_${RANDSTRING}.txt"

cat ${TFCSV} | awk -F $',' 'BEGIN {OFS = FS} $4 == "Yes" {print}' | cut -d"," -f1 > ${TMP_TF_FILE}
gunzip -c ${GTFFILE} > ${TMP_GTF_UNGZ}
grep -Fwf ${TMP_TF_FILE} ${TMP_GTF_UNGZ} | awk -F$"\t" 'BEGIN {OFS = FS} $3 == "gene" {print $1, $4, $5, $9}' > ${TMP_SEL_FILE}

echo -e 'ensembl_id\tchrom\tstart\tend\tname' > $OUTFILE
awk -f reformat_gtf_selection.awk ${TMP_SEL_FILE} >> $OUTFILE

rm -f ${TMP_TF_FILE}
rm -f ${TMP_GTF_UNGZ}
rm -f ${TMP_SEL_FILE}
