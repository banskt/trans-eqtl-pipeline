#!/bin/bash


PYENV="/usr/users/fsimone/opt/miniconda/3/envs/env3.6/bin/python"
methods="tejaas_perm tejaas_rand_perm matrixeqtl matrixeqtl_rand"
sbs="0.01 0.05"


for method in $methods; do
    if [ "${method}" == "tejaas_perm" ] || [ "${method}" == "tejaas_rand_perm" ]; then
        for sb in $sbs; do
            echo bsub -n 8 -R cbscratch -q mpi-long+ -x -W "240:00" \
                    -o "${method}_${sb}.out" \
                    -e "${method}_${sb}.err" \
                    $PYENV validate.py --method $method --sb $sb
        done;
    else
        echo bsub -n 8 -R cbscratch -q mpi-long+ -x -W "240:00" \
                -o "${method}.out" \
                -e "${method}.err" \
                $PYENV validate.py --method $method --sb 00
    fi
done;