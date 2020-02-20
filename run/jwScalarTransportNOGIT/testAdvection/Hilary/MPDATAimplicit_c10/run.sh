#!/bin/bash -e

rm -rf constant/polyMesh [0-9]*

blockMesh
mkdir 0
cp constant/Tsave constant/T_analytic_init 
setAnalyticTracerField
setVelocityField

cp 0/T_analytic 0/T

MPDATAadvectionFoam >& log & sleep 0.01; tail -f log

time=600
gmtFoam -time $time T
gv $time/T.pdf &

gmtFoam vT -time $time
gv $time/vT.pdf &


