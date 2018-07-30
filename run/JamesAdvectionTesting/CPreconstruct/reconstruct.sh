#!/bin/bash
set -e
rm -rf [0-9]*
blockMesh
mkdir 0
setInitialTracerField
testCPreconstruction
gmtFoam -time 0,1 mesh
