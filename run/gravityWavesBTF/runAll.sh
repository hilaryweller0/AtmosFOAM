#!/bin/bash -e

# create and run the gravity waves over orography test case from Schar

# mesh generation
blockMesh
add2dMountain
plotPatchData meshZoom
gv constant/meshZoom.eps &
createSpongeLayer

# initial conditions
mkdir -p 0
cp init_0/Exner_init constant
cp init_0/theta_init constant
cp init_0/Uf_init 0/Uf

# set temperature profile
setTheta

# set initialisation for H version
setExnerBalancedH

rm constant/*init

# run case
exnerFoamH |& tee log

# plot results
plotPatchData -time 18000 w
gv 18000/w.pdf &



