#!/bin/bash
set -e
rm -rf [0-9]*
blockMesh
mkdir 0
#cp init_0/Tf 0/Tf
setInitialTracerField
cp init_0/Uf 0/Uf
advectionConservativeF
gmtFoam -time 0,1 Tf
gmtFoam -time 1 TfRaw
zeroVerticalFaces Tf
cp 0/Tf 0/TfRaw
globalSum Tf
