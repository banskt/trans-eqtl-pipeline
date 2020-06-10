#!/bin/bash
#SBATCH -p hh
#SBATCH -t 0-12:00:00
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -J sboptim_msg
#SBATCH -o sboptim_msg.out
#SBATCH -e sboptim_msg.err

module load intel/compiler/64/2017/17.0.2
module load intel/mkl/64/2017/2.174
export OMP_NUM_THREADS=8
MPYTHON=/usr/users/fsimone/opt/miniconda/3/envs/env3.6/bin/python

RUN_PATH=/usr/users/fsimone/saikat/sb_optimization/sigbeta_optimizer.py
GXPRFILE=/cbscratch/franco/trans-eqtl/new_preprocess_feb2020_freeze/gtex_v8/expression/tpms/msg_tpms_qcfilter.txt.protein_coding_lncRNA_filtered
OUTFILE="/usr/users/fsimone/saikat/sb_optimization/gtex_v8/msg_sboptim.txt"
SAMPLE="/cbscratch/franco/datasets/gtex_v8/genotypes/gtex_v8_eur.sample"

${MPYTHON} ${RUN_PATH} --gx ${GXPRFILE} --out ${OUTFILE} --nsnp 5000 --samples ${SAMPLE}
