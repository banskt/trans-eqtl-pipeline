#!/bin/bash

FILE=$1
for variable in `grep -v -e '^#\|^if\|^fi' ${FILE} | grep -e "=" | sed '/^\s*$/d' | cut -d"=" -f1`; do
    unset $variable;
done
