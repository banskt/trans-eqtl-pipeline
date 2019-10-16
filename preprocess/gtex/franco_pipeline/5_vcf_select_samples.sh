#!/bin/bash

source "./CONFIG.sh"

while IFS='' read -r line || [ -n "$line" ]; do
    fullname=$(echo "$line" | cut -f 1 )
    shortname=$(echo "$line" | cut -f 2 )
    base=$(echo $fullname | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

    if [ $shortname = "wb" ] || [ $shortname = "hlv" ] || [ $shortname = "ms" ]  ; then
	    SAMPLEFILE=$EXPROUTDIR/$base.samples
	    VCFOUTDIR="$BASEDIR/genotypes/vcfs_$shortname"
	    mkdir -p $VCFOUTDIR
	    if [ -e $SAMPLEFILE ] ; then
	    	for chrom in `seq 1 22`; do
		    	INFILE="$VCFDIR/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_imput_info04_PASS_maf01_HWEp1E6_dbSNP135_ConstrVarIDs_chr${chrom}.vcf.gz"
		    	OUTFILE="$VCFOUTDIR/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_imput_info04_PASS_maf01_HWEp1E6_dbSNP135_ConstrVarIDs_${base}_"
		        $ENV $HOME/genomic_tools/select_samples_from_vcf.py --input ${INFILE} --outprefix ${OUTFILE} --incl-samples ${SAMPLEFILE}
		    done
	    else
	        echo "$fullname samplefile not found"
	    fi
	fi
done < $TISSUEFILE

