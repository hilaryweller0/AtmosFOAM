#!/bin/bash -e
#
# Create and run an OpenFOAM case for incompressible flow over orography

# clear out old stuff
rm -rf [0-9]* constant/polyMesh core log

# create mesh and plot
blockMesh
add2dMountain
gmtFoam mesh
gv constant/mesh.pdf &

# initialise the velocity field and the pressure
mkdir 0
cp init_0/* 0

# Solve incompressible Euler equations
icoFoamH >& log &
tail -f log

# plot the results
time=100
gmtFoam pU -time $time
gv $time/pU.pdf &

