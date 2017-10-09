#!/bin/bash
set -e
rm -rf [0-9]* processor*
blockMesh
terrainFollowingMesh

setTheta -CP

setExnerBalancedH -noInterpolate

cp init_0/Uf 0/Uf
createSpongeLayer

decomposePar -force -constant
for CASE in processor*; do
	fixProcessorFaceVelocities -case $CASE
done
set +e
mpirun --hostfile machines -np 2 exnerFoamCPinterpGrad -parallel
./diff.sh
