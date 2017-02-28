#!/bin/bash
set -e
rm -rf processor*
blockMesh
gmtFoam -time 0 mesh
gv 0/mesh.pdf &
echo SERIAL
testExtendedCentredFaceToFaceStencil
echo PARALLEL
decomposePar -force
mpirun -np 2 testExtendedCentredFaceToFaceStencil -parallel
mkdir -p processor0/constant/gmtDicts
mkdir -p processor1/constant/gmtDicts
cp constant/gmtDicts/mesh processor0/constant/gmtDicts/
cp constant/gmtDicts/mesh processor1/constant/gmtDicts/
gmtFoam -case processor0 -time 0 mesh
gmtFoam -case processor1 -time 0 mesh
gv processor0/0/mesh.pdf &
gv processor1/0/mesh.pdf &
