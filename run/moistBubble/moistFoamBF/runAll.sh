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
initMoistFoamBF
ln -sf ../0/Exner constant/ExnerRef

time=0
gmtFoam -time $time Exner
gv $time/Exner.pdf &
gmtFoam -time $time theta
gv $time/theta.pdf &
gmtFoam -time $time thetae
gv $time/thetae.pdf &
gmtFoam -time $time rv
gv $time/rv.pdf &
gmtFoam -time $time rl
gv $time/rl.pdf &

gmtFoam -time $time condenseRate
gv $time/condenseRate.pdf &

# run moistFoam
moistFoamBF >& log &
tail -f log

# debugging (differences from initial conditions
time=10
for var in thetae rv rl Exner theta; do
    echo $var
    sumFields $time ${var}Diff $time $var 0 ${var} -scale1 -1 | grep Min
    gmtFoam -time $time ${var}Diff > /dev/null
    gv $time/${var}Diff.pdf &
done
#gmtFoam -time $time condenseRate


# contours of vertical velocity
writeuvw Uf -time $time
gmtFoam -time $time w_BF02
gv $time/w_BF02.pdf &

writeuvw Uf -time $time
gmtFoam thetae -time $time
gv $time/thetae.pdf &

writeuvw Uf
gmtFoam thetae
eps2gif thetae.gif 0/thetae.pdf ??0/thetae.pdf 1000/thetae.pdf

