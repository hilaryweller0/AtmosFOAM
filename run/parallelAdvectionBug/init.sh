#!/bin/bash
set -e
rm -rf [0-9]* constant/polyMesh processor* *.dat *.cpt constant/*Terms
blockMesh
setInitialTracerField
setVelocityField
decomposePar
mpirun --hostfile machines -np 2 advectionFoam -parallel
reconstructPar
cp 0/T 0/T.par
mv 1/T 1/T.par
mv 2/T 2/T.par
advectionFoam
sumFields -scale0 1 -scale1 -1 1 T_diff 1 T 1 T.par
sumFields -scale0 1 -scale1 -1 2 T_diff 2 T 2 T.par
gmtFoam -time 1 T_diff
gv 1/T_diff.pdf &

