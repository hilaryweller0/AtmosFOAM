#!/bin/bash -e

# clear out old stuff
rm -r [0-9]* constant/polyMesh core log

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
gmtFoam -time $time Exner
gv $time/Exner.pdf &
gmtFoam -time $time p
gv $time/p.pdf &
gmtFoam -time $time T
gv $time/T.pdf &
gmtFoam -time $time thetaRho
gv $time/thetaRho.pdf &
gmtFoam -time $time thetae
gv $time/thetae.pdf &
gmtFoam -time $time qv
gv $time/qv.pdf &
gmtFoam -time $time ql
gv $time/ql.pdf &
gmtFoam -time $time Lv
gv $time/Lv.pdf &
gmtFoam -time $time es
gv $time/es.pdf &

gmtFoam -time $time condenseRate
gv $time/condenseRate.pdf &

# run moistFoam
moistFoam >& log &
tail -f log


# Differences from initial conditions
time=10
sumFields $time UfDiff $time Uf 0 Uf -scale1 -1 | grep Min
for var in Exner qv ql thetaRho ; do
    echo $var
    #sumFields $time ${var}Diff $time $var constant ${var}Ref -scale1 -1 | grep Min
    sumFields $time ${var}Diff $time $var 0 ${var} -scale1 -1 | grep Min
    gmtFoam -time $time ${var}Diff > /dev/null
    gv $time/${var}Diff.pdf &
done

# contours of vertical velocity
writeuvw Uf -time $time
gmtFoam -time $time w_BF02
gv $time/w_BF02.pdf &

writeuvw Uf -time $time
gmtFoam thetae -time $time
gv $time/thetae.pdf &

writeuvw Uf
gmtFoam thetae
eps2gif thetae.gif 0/thetae.pdf ???/thetae.pdf 1000/thetae.pdf

