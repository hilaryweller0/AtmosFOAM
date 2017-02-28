#!/bin/bash
set -e
rm -rf processor*
blockMesh
gmtFoam -time 0 mesh
gv 0/mesh.pdf &
testExtendedCentredFaceToFaceStencil
