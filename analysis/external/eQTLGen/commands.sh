#!/bin/bash

INPUT="2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt"
./sort_eqtlgen.sh ${INPUT} eQTLgen_significant_trans-eqtls_sorted.txt
./unique_eqtlgen.sh eQTLgen_significant_trans-eqtls_sorted.txt eQTLgen_unique_trans-eqtls_sorted.txt
tail -n +2 eQTLgen_unique_trans-eqtls_sorted.txt | awk '{printf "chr%d\t%d\t%d\n",$2,$3,$3+1}' > eQTLgen_bed_hg19.txt
## Next uplift the eQTLgen_bed_hg19.txt from https://genome.ucsc.edu/cgi-bin/hgLiftOver and get eQTLgen_bed_hg38_hglft.txt
#./hg38_convert_unique_eqtlgen.sh eQTLgen_unique_trans-eqtls_sorted.txt eQTLgen_bed_hg38_hglft.txt eQTLgen_unique_trans-eqtls_sorted_hg38.txt
