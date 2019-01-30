#!/bin/bash

FILE=$1

while read LINE; do
    echo $LINE | cut -d" " -f10- | sed 's/\ /\n/g' | cut -d":" -f1 | sort | uniq -c | grep "\./\."
done < $FILE
