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
initMoistFoam
ln -sf ../0/Exner constant/ExnerRef

time=0
for var in cpml es Exner kappam Lv T p ql qv qvs rho Rm thetae thetaRho; do
    gmtFoam -time $time $var
    gv $time/$var.pdf &
done
gmtFoam -time $time condenseRate
gv $time/condenseRate.pdf &

# run moistFoam
moistFoam >& log &
tail -f log

# debugging (differences from initial conditions
time=2
for var in  Exner ql qv thetaRho; do
    echo $var
    sumFields $time ${var}Diff $time $var 0 ${var} -scale1 -1 | grep -i min
    gmtFoam -time $time ${var}Diff > /dev/null
    gv $time/${var}Diff.pdf &
done

# absolute fields
time=2
for var in heatSourse condenseRate condForce condLim; do
    gmtFoam -time $time $var
    gv $time/$var.pdf &
done


# contours of vertical velocity
writeuvw Uf -time $time
gmtFoam -time $time w_BF02
gv $time/w_BF02.pdf &

moistThermoVars -time $time
writeuvw Uf -time $time
gmtFoam thetae -time $time
gv $time/thetae.pdf &

writeuvw Uf
gmtFoam thetae
eps2gif thetae.gif ?/thetae.pdf ??/thetae.pdf ???/thetae.pdf 1000/thetae.pdf

