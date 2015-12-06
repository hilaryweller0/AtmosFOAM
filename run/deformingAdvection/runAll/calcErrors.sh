#!/bin/bash -ve

if [ "$#" -ne 1 ]
then
   echo usage: post.sh case
   exit
fi

case=$1

# Difference between the numerical and analytic solutions
for time in $case/[0-9]*; do
    t=`filename $time`
    sumFields -case $case $t Tdiff $t T ../analytic/$t T -scale0 1 -scale1 -1
done

# Calculate error metrics
globalSum -case $case T
globalSum -case $case Tdiff

echo '#time l1 l2 linf mean var min max' > $case/errorNorms.dat
paste $case/globalSumTdiff.dat $case/../analytic/globalSumT.dat \
      $case/globalSumT.dat | tail -n +2 |\
    awk '{print $1, $2/$10, $3/$11, $4/$12, $5/$13, $6/$14, $23, $24}' \
    >> $case/errorNorms.dat

