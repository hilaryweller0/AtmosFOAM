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
ev $time/vorticityMesh.pdf

# Calculate the height in balance and plot
setBalancedHeightRC
gmtFoam -time $time hUmesh
gv $time/hUmesh.pdf &

# create initial conditions for partitioned equations
mv 0/Uf 0/bouyant.Uf
cp init_0/Uf 0/stable.Uf
rm 0/U 0/hU 0/p 0/streamFunction 0/velPot
setFields
sumFields 0 stable.sigma init_0 stable.sigma 0 buoyant.sigma -scale1 -1
gmtFoam -time 0 sigma
gv 0/sigma.pdf &

# Solve the SWE
partitionedShallowWaterFoam >& log & sleep 0.01; tail -f log

# Plots
time=10000
for plot in hu sigma; do
    gmtFoam -time $time $plot
    gv $time/$plot.pdf &
done


time=100000
case=.
postProcess -func vorticity -time $time -case $case
writeuvw -time $time vorticity -case $case
mv $case/$time/vorticityz $case/$time/vorticity
rm $case/$time/vorticity?
gmtFoam -time $time vorticity -case $case
gv $case/$time/vorticity.pdf &

# Animattion of vorticity
case=.
postProcess -func vorticity -case $case
writeuvw vorticity -case $case
for time in [0-9]*; do
    mv $case/$time/vorticityz $case/$time/vorticity
    rm $case/$time/vorticity?
done
gmtFoam vorticity -case $case
eps2gif vorticity.gif 0/vorticity.pdf ??????/vorticity.pdf \
                                      ???????/vorticity.pdf

# Make links for animategraphics
field=vorticity
DT=100000
mkdir -p animategraphics
for time in [0-9]*; do
    let t=$time/$DT
    ln -s ../$time/$field.pdf animategraphics/${field}_$t.pdf
done
