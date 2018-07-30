#!/bin/bash -e

# Clean up
##########
rm -rf [0-9]* constant/polyMesh

# Mesh and sponge generation
####################
blockMesh
terrainFollowingMesh
createSpongeLayer

# initial conditions
####################
cp -r init_0 0
cp init_0/theta constant/theta_init
cp init_0/Exner constant/Exner_init

# set temperature profile and hydrostatically balanced pressure
setTheta
setExnerBalanced
sed -i 's/fixedValue;/fixedFluxBuoyantExner; gradient uniform 0;/g' 0/Exner

# Plot initial theta
#gmtFoam -time 0 theta
#gv 0/theta.pdf &

# Run the case
##############
exnerFoam >& log & sleep 0.01; tail -f log

# Plot results
##############

time=18000
gmtFoam -time $time thetaU
gv $time/thetaU.pdf &

writeuvw -time $time Uf
mv $time/Ufz $time/Uz
gmtFoam -time $time w
gv $time/w.pdf &

sumFields $time ExnerDiff $time Exner 0 Exner -scale0 1 -scale1 -1
sumFields $time thetaDiff $time theta 0 theta -scale0 1 -scale1 -1
sumFields $time UDiff $time Uf 0 Uf -scale0 1 -scale1 -1
gmtFoam -time $time thetaUdiff
gv $time/thetaUdiff.pdf &

