#!/bin/bash -e

# clear out old stuff
rm -rf [0-9]* constant/polyMesh core log

# create mesh
blockMesh

# hydrostatically balanced initial conditions
rm -rf [0-9]* core
mkdir 0
cp -r init_0/* 0
setExnerBalancedH

# change Exner BC from fixedValue to hydroStaticExner
sed -i 's/fixedFluxBuoyantExner/partitionedHydrostaticExner/g' 0/Exner

# Add a warm perturnation
cp 0/theta 0/theta_init
makeHotBubble

# Partition into stable and buoyant fluids
mv 0/theta 0/buoyant.theta
mv 0/theta_init 0/stable.theta
mv 0/Uf 0/stable.Uf
cp 0/stable.Uf 0/buoyant.Uf
rm 0/thetaf

# create initial conditions
setFields
sumFields 0 stable.sigma init_0 stable.sigma 0 buoyant.sigma -scale1 -1

# Plot initial conditions
gmtFoam theta
gv 0/theta.pdf &
gmtFoam sigma
gv 0/sigma.pdf &

# Solve Euler equations
partitionedExnerFoam >& log & sleep 0.01; tail -f log

# animate the results
gmtFoam theta
animate 0/theta.pdf ???/theta.pdf 1000/theta.pdf

# Differences between non-partitioned run
sumFields 2 ExnerDiff 2 Exner ../standard/2 Exner -scale1 -1
sumFields 2 thetaDiff 2 stable.theta ../standard/2 theta -scale1 -1
sumFields 2 UfDiff 2 stable.Uf ../standard/2 Uf -scale1 -1
gmtFoam -time 2 ExnerDiff
gv 2/ExnerDiff.pdf &
gmtFoam -time 2 thetaDiff
gv 2/thetaDiff.pdf &

sumFields 2 fluxDiff 2 fluxSum ../standard/2 U -scale1 -1
sumFields 2 rhoDiff 2 stable.sigmaRho ../standard/2 rho -scale1 -1
sumFields 2 uDiff 2 stable.u ../standard/2 u -scale1 -1

sumFields 2 buoyant.sigma 2 stable.sigma init_0 stable.sigma -scale0 -1
gmtFoam -time 2 sigma
gv 2/sigma.pdf &

