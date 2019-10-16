#!/bin/bash

source "./CONFIG.sh"
CONVERT="/usr/users/fsimone/Tejaas_Pipeline/preprocess/convert_vcf_to_dosage.py"

while IFS='' read -r line || [ -n "$line" ]; do
    fullname=$(echo "$line" | cut -f 1 )
    shortname=$(echo "$line" | cut -f 2 )
    base=$(echo $fullname | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

    if [ $shortname = "wb" ] || [ $shortname = "hlv" ] || [ $shortname = "ms" ] ; then
	    VCFTISSUEDIR="$BASEDIR/genotypes/vcfs_$shortname"
	    DOSAGEOUTDIR="$BASEDIR/genotypes/dosages_$shortname"
	    mkdir -p $DOSAGEOUTDIR
    	for chrom in `seq 1 22`; do
    		INFILE="$VCFTISSUEDIR/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_imput_info04_PASS_maf01_HWEp1E6_dbSNP135_ConstrVarIDs_${base}_chr${chrom}.vcf.gz"
    		if [ -e $INFILE ]; then
	    		OUTFILE="$DOSAGEOUTDIR/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_imput_info04_PASS_maf01_HWEp1E6_dbSNP135_ConstrVarIDs_${base}_dosages_chr${chrom}.gz"
	    		echo bsub -n 4 -q mpi -a openmp -R cbscratch \
					-R span[hosts=1] \
					-o ${DOSAGEOUTDIR}/${chrom}.log \
					-e ${DOSAGEOUTDIR}/${chrom}.err \
					$ENV $PYTHON ${CONVERT} --vcf ${INFILE} --out ${OUTFILE}
			fi
	    done
	fi
done < $TISSUEFILE

VCFALLDIR="$BASEDIR/genotypes/vcfs_allsamples"
DOSAGEOUTDIR="$BASEDIR/genotypes/dosages_allsamples"
mkdir -p $DOSAGEOUTDIR
for chrom in `seq 1 22`; do
	INFILE="$VCFALLDIR/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_imput_info04_PASS_maf01_HWEp1E6_dbSNP135_ConstrVarIDs_chr${chrom}.vcf.gz"
	if [ -e $INFILE ]; then
		OUTFILE="$DOSAGEOUTDIR/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_imput_info04_PASS_maf01_HWEp1E6_dbSNP135_ConstrVarIDs_dosages_chr${chrom}.gz"
		bsub -n 4 -q mpi -a openmp -R cbscratch \
			-R span[hosts=1] \
			-o ${DOSAGEOUTDIR}/${chrom}.log \
			-e ${DOSAGEOUTDIR}/${chrom}.err \
			$ENV $PYTHON ${CONVERT} --vcf ${INFILE} --out ${OUTFILE}
	fi
done