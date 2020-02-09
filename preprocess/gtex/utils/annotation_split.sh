#!/bin/bash
echo "Pre-processing lookup table..."
zcat ${SRCLOOKUP} > ${GTOUT_ALL}/tmp_lookup.txt
sed -i '1d' ${GTOUT_ALL}/tmp_lookup.txt
rm -rf ${GTOUT_ALL}/${GTFILE_BASENAME}_chr*.annot

# Split the lookup file into 22 fil2s
awk '{printf "%s %s\n", $1, $7 >> "tmp_"$2".annot"}' ${GTOUT_ALL}/tmp_lookup.txt

# Rename the files properly
for CHRM in {1..22}; do
    mv tmp_chr${CHRM}.annot ${GTOUT_ALL}/${GTFILE_BASENAME}_chr${CHRM}.annot
done

rm -rf ${GTOUT_ALL}/tmp_lookup.txt

echo "Done."
