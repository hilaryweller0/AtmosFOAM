#!/bin/bash
set -e
rm -rf [0-9]*
blockMesh
setTheta -CP
setExnerBalancedH -noInterpolate
createSpongeLayer
#thermoVars
cp constant/Uf 0/Uf
exnerFoamCP
