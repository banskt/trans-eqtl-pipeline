#!/bin/bash

source "./CONFIG.sh"

while IFS='' read -r line || [ -n "$line" ]; do
    fullname=$(echo "$line" | cut -f 1 )
    shortname=$(echo "$line" | cut -f 2 )
    base=$(echo $fullname | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

    EXPRFILE="${EXPROUTDIR}/gtex.normalized.expression.${shortname}.txt"
    if [ -e $EXPRFILE ] ; then
        cut -f 2- ${EXPRFILE} |head -n 1 |tr "\t" "\n" > $EXPROUTDIR/$base.samples
    else
        echo "$fullname expression not found"
    fi
done < $TISSUEFILE

