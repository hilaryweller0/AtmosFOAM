#!/bin/bash
set -e
rm -r [0-9]*
blockMesh
setInitialTracerField
setVelocityField

advectionFoam

setAnalyticTracerField
sumFields -scale0 1 -scale1 -1 10000 T_diff 10000 T 10000 T_analytic

globalSum -time 10000 T_diff
mv globalSumT_diff.dat 10000
tail -n1 10000/globalSumT_diff.dat | cut -d' ' -f3 > 10000/l2errorT_diff.txt

globalSum -time 10000 T_analytic
mv globalSumT_analytic.dat 10000
tail -n1 10000/globalSumT_analytic.dat | cut -d' ' -f3 > 10000/l2errorT_analytic.txt

python3 -c "print(`paste -d'/' 10000/l2errorT_diff.txt 10000/l2errorT_analytic.txt`)" > 10000/l2errorT.txt
