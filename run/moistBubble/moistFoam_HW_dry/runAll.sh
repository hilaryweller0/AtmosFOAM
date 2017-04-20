#!/bin/bash -e

# clear out old stuff
rm -rf [0-9]*[0-9] constant/polyMesh core log

# create mesh
blockMesh
#gmtFoam mesh
#gv constant/mesh.pdf &

# create initial conditions
rm -rf [0-9]* core
mkdir 0
cp -r init_0/* 0
initMoistFoam_HW

sed -i 's/fixedValue;/hydrostaticExner; gradient uniform 0;/g' 0/Exner

# run moistFoam
moistFoam_HW >& log &
tail -f log

time=0
for var in airVapourRho waterVapourRho waterLiquidFrac theta Exner
do
    gmtFoam -time $time $var
    gv $time/$var.pdf &
done

# de-bugging
time=2
for var in theta airVapourRho waterLiquidFrac waterVapourRho; do
    sumFields $time ${var}Diff $time $var 0 ${var} -scale1 -1
    gmtFoam -time $time ${var}Diff
    gv $time/${var}Diff.pdf &
done
#for var in divu condensationOfwater heatSource; do
for var in theta airVapourRho waterLiquidFrac waterVapourRho; do
    gmtFoam -time $time $var
    gv $time/$var.pdf &
done

# output
writeuvw Uf -time $time
moistThermoVars_HW -time $time
gmtFoam thetae -time $time
gv $time/thetae.pdf &

#animation
writeuvw Uf
moistThermoVars_HW
gmtFoam thetae
eps2gif thetae.gif 0/thetae.pdf ???/thetae.pdf 1000/thetae.pdf


