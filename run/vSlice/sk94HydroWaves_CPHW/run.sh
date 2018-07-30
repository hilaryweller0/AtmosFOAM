#!/bin/bash -e
rm -rf [0-9]* constant/polyMesh/* log* *.dat
mkdir -p constant/polyMesh
ln -sf ../blockMeshDict constant/polyMesh/blockMeshDict
blockMesh
setTheta -CP
cp 0/theta constant/thetaRef
cp 0/thetaf constant/thetafRef
setExnerBalancedH -noInterpolate
perturbField
cp constant/Uf_init 0/Uf
exnerFoamCP >& log &
tail -f log

