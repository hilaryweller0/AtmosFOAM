#!/bin/bash -e

# clear out old stuff
rm -rf [0-9]* constant/polyMesh core log

# create mesh and plot
blockMesh
gmtFoam mesh
gv constant/mesh.pdf &

# hydrostatically balanced initial conditions
rm -rf [0-9]* core
mkdir 0
cp -r init_0/* 0
setExnerBalancedH

# change Exner BC from fixedValue to hydroStaticExner
sed -i 's/fixedValue;/fixedFluxBuoyantExner; gradient uniform 0;/g' 0/Exner

# Add a warm perturnation
cp 0/theta 0/theta_init
makeHotBubble

# Plot initial potential temperature
gmtFoam theta
gv 0/theta.pdf &

# Solve Euler equations
exnerFoamH >& log &
tail -f log

# animate the results
gmtFoam theta
animate 0/theta.pdf ???/theta.pdf 1000/theta.pdf

