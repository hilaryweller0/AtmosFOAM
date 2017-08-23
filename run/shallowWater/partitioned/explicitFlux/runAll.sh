#!/bin/bash -e

# Create the case, run and post-process

# clear out old stuff
rm -rf [0-9]*[0-9] constant/polyMesh core log legends gmt.history

# Create mesh
blockMesh

# Create initial conditions
rm -rf 0 [0-9]*[0-9] core
cp -r init_0 0
time=0
# Create Gaussian patches of voriticty
setGaussians initDict
# Invert to find the wind field
invertVorticity -time $time initDict

gmtFoam -time $time vorticityMesh
#ev $time/vorticityMesh.pdf

# Calculate the height in balance and plot
setBalancedHeightRC
gmtFoam -time $time hUmesh
#gv $time/hUmesh.pdf &

# create initial conditions for partitioned equations
mv 0/Uf 0/buoyant.Uf
cp init_0/Uf 0/stable.Uf
rm 0/U 0/hU 0/p 0/streamFunction 0/velPot
setFields
sumFields 0 stable.sigma init_0 stable.sigma 0 buoyant.sigma -scale1 -1
gmtFoam -time 0 sigma
#gv 0/sigma.pdf &

# Solve the SWE
partitionedShallowWaterFoamFluxExp >& log & sleep 0.01; tail -f log

# Plots
time=10000
for plot in hu sigma; do
    gmtFoam -time $time $plot
    gv $time/$plot.pdf &
done

gmtPlot plots/plotEnergy.gmt

# Plot sigma and hu for all times
for plot in hu sigma; do
    gmtFoam $plot
    eps2gif $plot.gif 0/$plot.pdf ?????/$plot.pdf ??????/$plot.pdf
done

