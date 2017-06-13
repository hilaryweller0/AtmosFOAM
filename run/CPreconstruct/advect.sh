#!/bin/bash
set -e
rm -rf [0-9]*
blockMesh
mkdir 0
cp init_0/Tf 0/Tf
cp init_0/Uf 0/Uf
advectionConservativeF
zeroVerticalFaces Tf
globalSum Tf
gmtFoam -time 0,1 mesh
