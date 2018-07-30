#!/bin/bash
set -e
rm -rf [0-9]*
blockMesh
setTheta -CP
setExnerBalancedH -noInterpolate
createSpongeLayer
setInitialTracerField
mv 0/T constant/radiation
mv 0/Tf constant/radiationf
gmtFoam -time constant radiation
cp constant/Uf 0/Uf
exnerFoamCPinterpGrad
