#!/bin/bash -e
#
# Create and run an incompressible flow over backward facing step using
# icoFoamH

# clear out old stuff
rm -rf [0-9]* constant/polyMesh core log

# create mesh and plot
blockMesh
gmtFoam mesh
gv constant/mesh.pdf &

# initialise the velocity field and the pressure
mkdir 0
cp init_0/* 0

# Solve incompressible Euler equations
icoFoamH >& log &
tail -f log

# plot the results
time=0.001
gmtFoam pU -time $time
gv $time/pU.pdf &


gmtFoam pU
eps2gif pU.gif 0.??/pU.pdf &

