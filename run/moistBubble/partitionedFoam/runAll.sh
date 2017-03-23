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
ss=0.9
bs=0.1
for time in 2; do
    for var in rho Exner theta Uf ; do
        sumFields $time ${var}Diff $time $var ../moistFoam_HW/$time $var -scale1 -1
        read -p "Press enter to continue"
    done
    sumFields $time fluxDiff $time flux ../moistFoam_HW/$time U -scale1 -1
    read -p "Press enter to continue"
    sumFields $time fluxDiff $time stable.flux ../moistFoam_HW/$time U -scale1 -$ss
    read -p "Press enter to continue"
    sumFields $time fluxDiff $time buoyant.flux ../moistFoam_HW/$time U -scale1 -$bs
done

# Differences between partitions
ss=0.9
bs=0.1
for time in 2; do
    for var in waterLiquidFrac waterVapourRho airVapourRho theta T Uf; do
        sumFields $time ${var}Diff $time stable.$var $time buoyant.$var -scale1 -1
        read -p "Press enter to continue"
    done
    for var in sigma flux; do
        sumFields $time ${var}Diff $time stable.$var $time buoyant.$var -scale0 $bs -scale1 -$ss
        read -p "Press enter to continue"
    done
    read -p "Press enter to continue to next time-step"
done

# Differences from initial conditions
time=10
sumFields $time thetaDiff $time theta 0 theta -scale1 -1
sumFields $time thetaDiff $time stable.theta 0 theta -scale1 -1
sumFields $time thetaDiff $time buoyant.theta 0 theta -scale1 -1
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


