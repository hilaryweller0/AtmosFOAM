#!/bin/bash -e

# clear out old stuff
rm -rf [0-9]* constant/polyMesh core log

# create mesh
blockMesh
#gmtFoam mesh
#gv constant/mesh.pdf &

# create initial conditions
rm -rf [0-9]* core
mkdir 0
cp -r init_0/* 0
initMoistFoam_HW

# run moistFoam
moistFoam_HW


