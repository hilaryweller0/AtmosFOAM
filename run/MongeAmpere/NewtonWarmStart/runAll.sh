#!/bin/bash -e

rm -rf core constant [0-9]*
blockMesh
mkdir -p constant/rMesh 0
cp init0/Phi 0
cp -r constant/polyMesh constant/rMesh
NEWTON2D-pab >& log &
tail -f log

# next edit system/controlDict so that endTime is 100 and it will restart from time step 50 
