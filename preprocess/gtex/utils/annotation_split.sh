#!/bin/bash
echo "Pre-processing lookup table..."
zcat ${SRCLOOKUP} > ${GTOUT_ALL}/tmp_lookup.txt
sed -i '1d' ${GTOUT_ALL}/tmp_lookup.txt
rm -rf ${GTOUT_ALL}/${GTFILE_BASENAME}_chr*.annot

# Split the lookup file into 22 files
awk '{printf "%s %s\n", $3, $6 >> "tmp_"$1".annot"}' ${GTOUT_ALL}/tmp_lookup.txt

# Rename the files properly
for CHRM in {1..22}; do
    mv tmp_${CHRM}.annot ${GTOUT_ALL}/${GTFILE_BASENAME}_chr${CHRM}.annot
done

rm -rf ${GTOUT_ALL}/tmp_lookup.txt

echo "Done."
