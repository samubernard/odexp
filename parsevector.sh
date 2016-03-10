#!/bin/bash

# x0:h:x1

typ=$1
file=$2

v=(`awk -F ' ' -v vartype=$typ '$1 ~ vartype && $2 ~ /\[[0-9]+\]/ {split($2,a,/[\[\]]/); print a[1]}' $file`)
nv=(`awk -F ' ' -v vartype=$typ '$1 ~ vartype && $2 ~ /\[[0-9]+\]/ {split($2,a,/[\[\]]/); print a[2]}' $file`)
ev=(`awk -F ' ' -v vartype=$typ '$1 ~ vartype && $2 ~ /\[[0-9]+\]/ { print $3}' $file`)

n=`expr ${#v[@]} - 1`
if [ "$n" -ge 0 ]
then
    for k in `seq 0 $n`
    do
        echo "    for (_i=0; _i<${nv[$k]}; _i++)" >>.odexp/model.c
        echo "    {" >>.odexp/model.c
        echo "      ${v[$k]}[_i] = ${ev[k]};" >>.odexp/model.c
        echo "    }" >>.odexp/model.c
    done
fi
