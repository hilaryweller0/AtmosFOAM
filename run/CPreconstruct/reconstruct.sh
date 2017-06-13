#!/bin/bash
set -e
rm -rf [0-9]*
blockMesh
mkdir 0
cp init_0/Tf 0/Tf
testCPreconstruction
gmtFoam -time 0,1 mesh
