#!/bin/bash
set -e
rm -rf 1 constant/polyMesh processor* *.dat *.cpt constant/*Terms
blockMesh
#setInitialTracerField
setVelocityField
decomposePar
mpirun --hostfile machines -np 2 advectionFoam -parallel -forwardEuler
reconstructPar
cp 0/T 0/T.par
mv 1/T 1/T.par
advectionFoam -forwardEuler
sumFields -scale0 1 -scale1 -1 1 T_diff 1 T 1 T.par
gmtFoam -time 1 T_diff
gv 1/T_diff.pdf &

