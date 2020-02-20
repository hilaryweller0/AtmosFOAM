#!/bin/bash -e

if [ "$#" -ne 2]
then
   echo usage: plotOne.sh dirName $time
   exit
fi

dir=$1
time=$2

rm -f $dir/$time/CoImp
if [[ $dir == *"implicitWhereNeeded"* ]]; then
    ln -s ../0/CoImp $dir/$time
fi

gmtFoam -case $dir -time $time T
if [[ $dir == *"MPDATA"* ]]; then
    gmtFoam -case $dir -time $time vT
fi

