#!/bin/bash -e

# Clear out old mesh and results
rm -rf constant/polyMesh [0-9]*

# Create mesh and initial conditions
blockMesh
mkdir 0
cp init_0/Tsave constant/T_analytic_init 
setAnalyticTracerField
setVelocityField
cp 0/T_analytic 0/T

# Plot initial conditions
gmtFoam -time 0 UT
gv 0/UT.pdf &

# Solve linear advection equation
scalarTransportFoamCN >& log & sleep 0.01; tail -f log

# Plot final results
time=600
gmtFoam -time $time T
gv $time/T.pdf &

