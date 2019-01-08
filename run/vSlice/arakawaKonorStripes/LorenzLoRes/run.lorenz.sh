#!/bin/bash
set -e
rm -rf [0-9]* processor*
blockMesh

setInitialTracerField
mv 0/T 0/theta.anom

#terrainFollowingMesh
setTheta

setExnerBalancedH

mv 0/theta 0/theta.bg

sumFields -scale0 1 -scale1 1 0 theta 0 theta.bg 0 theta.anom

cp init_0/Uf 0/Uf
createSpongeLayer
set +e

decomposePar -force -constant
mpirun --hostfile machines -np 2 exnerFoamH -parallel
./diff.sh