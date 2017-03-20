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

# run partitionedMoistFoam
partitionedMoistFoam >& log &
tail -f log

# Differences from non-partitioned code
time=2
sumFields $time Udiff $time Uf ../moistFoam_HW/$time Uf -scale1 -1
sumFields $time ExnerDiff $time Exner ../moistFoam_HW/$time Exner -scale1 -1

sumFields $time rhoDiff $time rho ../moistFoam_HW/$time rho -scale1 -1
sumFields $time thetaDiff $time stable.theta ../moistFoam_HW/$time theta -scale1 -1
sumFields $time waterVapourRhoDiff $time stable.waterVapourRho ../moistFoam_HW/$time waterVapourRho -scale1 -1
sumFields $time sigmaTmp $time stable.sigma $time stable.sigma -scale1 0

# Differences from initial conditions
time=10
sumFields $time thetaDiff $time theta 0 stable.theta -scale1 -1
sumFields $time thetaDiff $time stable.theta 0 stable.theta -scale1 -1
sumFields $time thetaDiff $time buoyant.theta 0 buoyant.theta -scale1 -1
gmtFoam -time $time thetaDiff
gv $time/thetaDiff.pdf &
sumFields $time ExnerDiff $time Exner 0 Exner -scale1 -1
gmtFoam -time $time ExnerDiff
gv $time/ExnerDiff.pdf &


time=100
gmtFoam -time $time thetaU
gv $time/thetaU.pdf &

# plotting output
writeuvw Uf -time $time
partitionedMoistThermoVars -time $time
gmtFoam thetae -time $time
gv $time/thetae.pdf &

#animation
writeuvw Uf
moistThermoVars_HW
gmtFoam thetae
eps2gif thetae.gif 0/thetae.pdf ???/thetae.pdf 1000/thetae.pdf


