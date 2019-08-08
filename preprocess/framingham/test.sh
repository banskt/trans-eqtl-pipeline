#!/bin/bash

COMMAND='cat /cbscratch/franco/datasets/FHS/genotypes/c1/chr22_c1.fhs.dosages.c1.sample <(tail -n +3 /cbscratch/franco/datasets/FHS/genotypes/c2/chr22_c2.fhs.dosages.c2.sample ) > /cbscratch/franco/datasets/FHS/genotypes/merged/chr22.fhs.dosages.sample'

echo $COMMAND
eval $COMMAND
