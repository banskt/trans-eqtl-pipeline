
TISSUES="aa as at esom esomu fib lu ms nt sse thy wb"
DATADIR="/cbscratch/franco/datasets"
GENEINFOFILE="${DATADIR}/GENCODE/gencode.v19.annotation.gtf.gz"
SAMPLE="/cbscratch/franco/datasets/gtex/gtex.sample"
for T in $TISSUES; do
    EXPR=/cbscratch/franco/trans-eqtl/new_preprocess/gtex/expression/lmcorrected/${T}_age_lmcorrected.precorr.txt
    python filter_gencode_expr.py --gx ${EXPR} \
                                  --donors ${SAMPLE} --dataset gtex \
                                  --gtf ${GENEINFOFILE} \
                                  --biotype protein_coding;
done;

