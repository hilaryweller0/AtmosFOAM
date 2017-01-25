#!/bin/bash
set -e
rm -rf [0-9]*
blockMesh
setTheta
setExnerBalancedH
createSpongeLayer
thermoVars
cp constant/Uf 0/Uf
exnerFoamH
