#!/bin/bash -e

# create the case and run the Newton solver for the MA equation

# clear out old results and old meshes
rm -rf [0-9]* constant/polyMesh constant/rMesh log

# create the computational mesh and the physical mesh (the rMesh)
blockMesh
mkdir constant/rMesh
cp -r constant/polyMesh constant/rMesh

# create the initial conditions
mkdir 0
cp constant/Phi_0 0/Phi

# Solve the MA equation
touch log
tail -f log &
Newton_HW > log 2>&1
pkill tail

# plots of Phi and of the new mesh
gmtFoam -latestTime Phi
evince */Phi.pdf &
gmtFoam -latestTime -region rMesh mesh
evince */mesh.pdf &

exit

time=100
gmtFoam -time $time Phi
evince $time/Phi.pdf &
gmtFoam -time $time -region rMesh mesh
evince $time/mesh.pdf &

