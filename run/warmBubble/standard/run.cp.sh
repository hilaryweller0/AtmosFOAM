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
setExnerBalancedH -noInterpolate

# change Exner BC from fixedValue to hydroStaticExner
sed -i 's/fixedValue;/fixedFluxBuoyantExner; gradient uniform 0;/g' 0/Exner

# Add a warm perturnation
cp 0/theta 0/theta_init
cp 0/thetaf 0/thetaf_init
makeHotBubble

# Plot initial potential temperature
gmtFoam thetaf
gv 0/thetaf.pdf &

# Solve Euler equations
exnerFoamCP >& log &
tail -f log

# animate the results
gmtFoam thetaf
animate 0/thetaf.pdf ???/thetaf.pdf 1000/thetaf.pdf

