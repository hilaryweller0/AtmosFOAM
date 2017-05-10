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
partitionedMoistFoam >& log & sleep 0.01; tail -f log

# Differences from non-partitioned code
ss=0.9
bs=0.1
for time in 2; do
    for var in Exner dRhodt Psi theta Uf ; do
        sumFields $time ${var}Diff $time $var ../moistFoam_HW/$time $var -scale1 -1
        read -p "That was $var, press enter to continue"
    done
    sumFields $time fluxDiff $time flux ../moistFoam_HW/$time U -scale1 -1
    read -p "That was the flux, Press enter to continue"
    sumFields $time fluxDiff $time stable.flux ../moistFoam_HW/$time U -scale1 -$ss
    read -p "Press enter to continue"
    sumFields $time fluxDiff $time buoyant.flux ../moistFoam_HW/$time U -scale1 -$bs
    
    for var in airVapourRho theta Uf waterLiquidFrac waterVapourRho; do
        sumFields $time ${var}Diff $time buoyant.$var ../moistFoam_HW/$time $var -scale1 -1
        read -p "That was $var for partition buoyant, press enter to continue"
    done

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
sumFields $time thetaDiff $time theta constant thetaRef -scale1 -1
sumFields $time thetaDiff $time stable.theta constant thetaRef -scale1 -1
sumFields $time thetaDiff $time buoyant.theta constant thetaRef -scale1 -1
gmtFoam -time $time thetaDiff
gv $time/thetaDiff.pdf &
sumFields $time ExnerDiff $time Exner 0 Exner -scale1 -1
gmtFoam -time $time ExnerDiff
gv $time/ExnerDiff.pdf &
sumFields $time waterLiquidFracDiff $time stable.waterLiquidFrac \
          0 stable.waterLiquidFrac -scale1 -1
gmtFoam -time $time waterLiquidFracDiff
gv $time/waterLiquidFracDiff.pdf &
sumFields $time sigmaDiff $time buoyant.sigma 0 buoyant.sigma -scale1 -1
gmtFoam -time $time sigmaDiff
gv $time/sigmaDiff.pdf &
gmtFoam -time $time rhoDiff
gv $time/rhoDiff.pdf &

# Different stable and buoyant partitions
for time in 0 ?? [1-3]??; do
    sumFields $time thetaDiff $time stable.theta constant thetaRef -scale1 -1
    gmtFoam -time $time stableThetaDiff
    sumFields $time thetaDiff $time buoyant.theta constant thetaRef -scale1 -1
    gmtFoam -time $time buoyantThetaDiff
#    sumFields $time ExnerDiff $time Exner 0 Exner -scale1 -1
#    gmtFoam -time $time ExnerDiff
#    gmtFoam -time $time sigma
done
eps2gif buoyantThetaDiff.gif 0/buoyantThetaDiff.pdf \
    ??/buoyantThetaDiff.pdf ???/buoyantThetaDiff.pdf
eps2gif stableThetaDiff.gif 0/stableThetaDiff.pdf \
    ??/stableThetaDiff.pdf ???/stableThetaDiff.pdf
#eps2gif ExnerDiff.gif 0/ExnerDiff.pdf ??/ExnerDiff.pdf ???/ExnerDiff.pdf
#eps2gif sigma.gif 0/sigma.pdf ??/sigma.pdf ???/sigma.pdf

time=800
for var in sigma sigmaRho theta waterLiquidFrac waterVapourRho; do
    for part in buoyant stable; do
        sumFields $time varDiff $time $part.$var 0 $part.$var -scale1 -1
        gmtFoam -time $time varDiff
        echo $part $var
        gv $time/varDiff.pdf
    done
done

# Differences from other code
time=100

sumFields $time thetaDiff $time stable.theta ../moistFoam_HW/$time

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


