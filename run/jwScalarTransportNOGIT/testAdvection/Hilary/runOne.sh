##!/bin/bash -e

if [ "$#" -ne 1 ]
then
   echo usage: runOne.sh dirName
   exit
fi

dir=$1

rm -rf $dir/constant/polyMesh $dir/[0-9]*

blockMesh -case $dir
mkdir $dir/0
cp $dir/constant/Tsave $dir/constant/T_analytic_init 
setAnalyticTracerField -case $dir
setVelocityField -case $dir

cp $dir/0/T_analytic $dir/0/T

if [[ $dir == *"MPDATA"* ]]; then
    MPDATAadvectionFoam -case $dir |&  tee $dir/log
else
    implicitExplicitAdvectionFoam -case $dir |&  tee $dir/log
fi

