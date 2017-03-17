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

time=2
sumFields $time Udiff $time uTmp ../moistFoam_HW/$time uTmp -scale1 -1
sumFields $time ExnerDiff $time Exner ../moistFoam_HW/$time Exner -scale1 -1
sumFields $time divUdiff $time divU ../moistFoam_HW/$time divU -scale1 -1
sumFields $time PsiDiff $time Psi ../moistFoam_HW/$time Psi -scale1 -1
sumFields $time gradPcoeff2Diff $time gradPcoeff2 ../moistFoam_HW/$time gradPcoeff2 -scale1 -1
gmtFoam -time $time Udiff
gv $time/Udiff.pdf &

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


