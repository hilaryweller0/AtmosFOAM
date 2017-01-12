#!/bin/bash
set -e
rm -rf 0/Tf 1 constant/polyMesh processor* *.dat *.cpt constant/*Terms
blockMesh
#setInitialTracerField
setVelocityField
decomposePar
mpirun --hostfile machines -np 2 advectionFoam -parallel -forwardEuler
#mpirun --hostfile machines -np 2 interpolateT -parallel
#mpirun --hostfile machines -np 2 xterm -e "gdb -ex run --args interpolateT -parallel"
reconstructPar
cp 0/T 0/T.par
mv 1/T 1/T.par
#mv 0/Tf 0/Tf.par
advectionFoam -forwardEuler
#interpolateT
sumFields -scale0 1 -scale1 -1 1 T_diff 1 T 1 T.par
#sumFields -scale0 1 -scale1 -1 0 Tf_diff 0 Tf 0 Tf.par
gmtFoam -time 1 T_diff
#gmtFoam -time 0 Tf_diff
xdg-open 1/T_diff.pdf &>/dev/null
#gv 0/Tf_diff.pdf &
