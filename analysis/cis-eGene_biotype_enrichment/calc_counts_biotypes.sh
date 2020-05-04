#!/bin/bash

INFILE="$1"
TFFILE="$2"

mapfile -t BIOTYPES < <( cut -f4 ${INFILE} | sort | uniq -c | sort -nr | awk '{print $2}' )
TISSUES=$( cut -f1 ${INFILE} | sort | uniq )
RAND=$( cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1 )

TFLIST="tflist_${RAND}_tmp.txt"
cut -f1 ${TFFILE} | tail -n +2 > ${TFLIST}

function result_string() {
  INFILE=$1
  TOUTFILE=$2
  TFLIST=$3
  TISSUE=$4
  RESSTRING="${TISSUE}\t"
  mapfile -t BIOTYPES < <( cut -f4 ${INFILE} | sort | uniq -c | sort -nr | awk '{print $2}' )
  for BIOTYPE in ${BIOTYPES[@]}; do
    BNUM="0"
    ISPRESENT=$( grep -cw "${BIOTYPE}" ${TOUTFILE})
    if [[ $ISPRESENT == 1 ]]; then
      BNUM=$( grep -w "${BIOTYPE}" ${TOUTFILE} | awk '{print $2}')
      #echo "$BIOTYPE $BNUM"
    fi
    RESSTRING+="${BNUM}\t"
  done
  TMP_TFCHECK="tmp_cisgenes_${TISSUE}_${RAND}.txt"
  if [ ! "${TISSUE}" == "all" ]; then
      awk -v MTISSUE=${TISSUE} '$1 == MTISSUE {print $3}' ${INFILE} | sort | uniq > ${TMP_TFCHECK}
  else
      cut -f3 ${INFILE} | sort | uniq > ${TMP_TFCHECK}
  fi
  NTF=$( grep -c -f ${TMP_TFCHECK} ${TFLIST} )
  rm -f ${TMP_TFCHECK}
  RESSTRING+="${NTF}"
  echo $RESSTRING
}

#HEADER=$( echo -e "TISSUE\t"; IFS=$'\t'; echo "${BIOTYPES[*]}" )
#HEADER+="\ttranscription_factors"
HEADER=$( echo "TISSUE" "${BIOTYPES[@]}" transcription_factors | tr " " "\t" )
echo -e ${HEADER}

for TISSUE in ${TISSUES}; do
  TOUTFILE=${TISSUE}_${RAND}_tmp.txt
  awk -v MTISSUE=${TISSUE} '$1 == MTISSUE {print $3, $4}' ${INFILE} | sort | uniq | cut -d" " -f2 | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > ${TOUTFILE}
  RESSTRING=$( result_string ${INFILE} ${TOUTFILE} ${TFLIST} ${TISSUE} )
  echo -e $RESSTRING
  rm -f ${TOUTFILE}
done

## Get the background
BGTMP="background_${RAND}_tmp.txt"
cut -f3,4 ${INFILE} | sort | uniq | awk '{print $2}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > ${BGTMP}
RESSTRING=$( result_string ${INFILE} ${BGTMP} ${TFLIST} all )
echo -e $RESSTRING
rm -f ${BGTMP}

rm -f ${TFLIST}
