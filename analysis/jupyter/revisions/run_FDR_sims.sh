#!/bin/bash -e

BASEDIR="/cbscratch/franco/trans-eqtl/simulation/gtex_v6_ldpruned_allTFstrength"
PYTHON37="/cbscratch/franco/envs/pyenv37/bin/python"

####### IMPORTANT ########
# parameters hardoced in the python script!
# - sigma_beta
# - parameters for each software and gene expression (raw, qn, cclm, etc)
# - algorithms to run


NTRANSS="50 100 150"
TF_TRANS_SHAPES="5 10 15 20"

for NTRANS in $NTRANSS; do
    for TF_TRANS_SHAPE in $TF_TRANS_SHAPES; do
        PARAMSTR="12639_450_10_800_30_${NTRANS}_100_0.01_0.5_0.0_1.0_0.6_4.0_0.1_${TF_TRANS_SHAPE}_0.02"

        PARAMDIR="${BASEDIR}/${PARAMSTR}" #12639_450_10_800_30_100_100_0.01_0.5_0.0_1.0_0.6_4.0_0.1_20_0.02
        OUTDIR="${PARAMDIR}/rocdata_franco_nocis"
        mkdir -p $OUTDIR

        STARTSIM=1
        NSIM=20
        ENDSIM=$(( STARTSIM + NSIM ))
        for (( SIM=$STARTSIM; SIM<$ENDSIM; SIM++ )); do

            SIMINDEX=`echo $SIM | awk '{printf "%03d", $1}'`
            SIMNAME="sim${SIMINDEX}"
            echo $SIMNAME
            # echo "${PYTHON37} calc_FDR_sims_eachSNP.py --input ${PARAMDIR} --simname ${SIMNAME} --outdir ${OUTDIR}"
            sbatch -p hh -N 1 -n 4 --mem=40G -t 1-00:00:00 -o "${OUTDIR}/${SIMNAME}.out" -e "${OUTDIR}/${SIMNAME}.err" \
                    --wrap="${PYTHON37} calc_FDR_sims_eachSNP.py --input ${PARAMDIR} --simname ${SIMNAME} --outdir ${OUTDIR}"
        done;
    done;
done;