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

cp init_0/U 0/U
createSpongeLayer
set +e

decomposePar -force -constant
#exnerFoamRC
mpirun -np 4 exnerFoamRC -parallel
./diff.sh
