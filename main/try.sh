#!/bin/bash

function foo() {
    local NAME=$1
    local SURNAME=$2
    echo "Hello $NAME $SURNAME"
}

function tejaas_fixpval() {
    local PYENV=$1
    local PYSCRIPT=$2
    local DIR=$3
    local CWD=$( pwd )
    echo $PYENV $PYSCRIPT $DIR
    #echo rr.txt rr_old.txt
    #${PYENV} ${PYSCRIPT} --infile $DIR/rr_old.txt --outfile $DIR/rr.txt
    cd $CWD
}

#foo Saikat Banerjee
#foo Rikhia Ghosh

tejaas_fixpval python some_script some_directory
