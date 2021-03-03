#!/bin/bash -e

BASEDIR="/cbscratch/franco/trans-eqtl/simulation/gtex_v6_ldpruned/12639_450_10_800_30_100_100_0.01_0.5_0.0_1.0_0.6_4.0_0.1_20_0.02"
PYTHON37="/usr/users/fsimone/opt/miniconda/3/envs/pyenv37/bin/python"
OUTDIR="${BASEDIR}/rocdata_franco_nocis"
mkdir -p $OUTDIR

STARTSIM=1
NSIM=20
ENDSIM=$(( STARTSIM + NSIM ))
for (( SIM=$STARTSIM; SIM<$ENDSIM; SIM++ )); do

    SIMINDEX=`echo $SIM | awk '{printf "%03d", $1}'`
    SIMNAME="sim${SIMINDEX}"
    echo $SIMNAME
    #echo "${PYTHON37} calc_tpr_fpr_sims.py --input ${BASEDIR} --simname ${SIMNAME} --outdir ${OUTDIR}"
    sbatch -p hh -N 1 -n 8 -t 2-00:00:00 -o "${OUTDIR}/${SIMNAME}.out" -e "${OUTDIR}/${SIMNAME}.err" \
            --wrap="${PYTHON37} calc_tpr_fpr_sims.py --input ${BASEDIR} --simname ${SIMNAME} --outdir ${OUTDIR}"
done;