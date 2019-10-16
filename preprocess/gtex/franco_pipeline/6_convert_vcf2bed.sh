#!/bin/bash

source "./CONFIG.sh"

while IFS='' read -r line || [ -n "$line" ]; do
    fullname=$(echo "$line" | cut -f 1 )
    shortname=$(echo "$line" | cut -f 2 )
    base=$(echo $fullname | sed 's/ - /_/g' | sed 's/ /_/g' | sed 's/(//g' | sed 's/)//g' )

    if [ $shortname = "wb" ] || [ $shortname = "hlv" ] || [ $shortname = "ms" ] ; then
        VCFINPUTDIR="$BASEDIR/genotypes/vcfs_split_$shortname"
        SAMPLEFILE="$EXPROUTDIR/$base.samples"
        if [ -e $SAMPLEFILE ] ; then
            BEDOUTDIR="$VCFINPUTDIR/bed"
            mkdir -p $BEDOUTDIR
            for chrom in `seq 1 22`; do
                INFILE="$VCFINPUTDIR/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_${base}_chr${chrom}.vcf.gz"
                OUTFILE="$BEDOUTDIR/GTEx_nomissing_chr${chrom}"
                $PLINK2 --vcf ${INFILE} dosage=DS --make-bed --maf 0.1 --max-maf 0.9 --geno dosage --hard-call-threshold 0.499999 --out ${OUTFILE} > ${OUTFILE}.log
            done
        else
            echo "$fullname samplefile not found"
        fi
    fi
done < $TISSUEFILE