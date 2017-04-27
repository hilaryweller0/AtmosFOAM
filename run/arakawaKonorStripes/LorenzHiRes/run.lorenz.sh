#!/bin/bash
set -e
rm -rf [0-9]* processor*
blockMesh

setInitialTracerField
mv 0/T 0/theta.anom

setTheta

setExnerBalancedH

mv 0/theta 0/theta.bg

sumFields -scale0 1 -scale1 1 0 theta 0 theta.bg 0 theta.anom

cp init_0/Uf 0/Uf
set +e
createSpongeLayer

decomposePar -force -constant
#mpirun --hostfile machines -np 2 exnerFoamH -parallel >& log &
#tail -f log
#reconstructPar

for t in [0-9]*
do
	sumFields -scale0 1 -scale1 -1 $t theta_diff $t theta 0 theta.bg
done

