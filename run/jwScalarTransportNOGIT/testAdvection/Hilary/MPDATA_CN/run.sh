#!/bin/bash -e

rm -rf constant/polyMesh [0-9]*

blockMesh
mkdir 0
cp constant/Tsave constant/T_analytic_init 
setAnalyticTracerField
setVelocityField

cp 0/T_analytic 0/T

#gmtFoam -time 0 UT
#gv 0/UT.pdf &

MPDATA_CN >& log & sleep 0.01; tail -f log

time=600
gmtFoam -time $time T
gv $time/T.pdf &

gmtFoam vT
eps2gif vT.gif 0/vT.pdf ??/vT.pdf ???/vT.pdf &

