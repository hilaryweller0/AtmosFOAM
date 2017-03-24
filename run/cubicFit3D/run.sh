#!/bin/bash
set -e
rm -rf [0-9]* processor* *dat *png
blockMesh
setInitialTracerField
setVelocityField
decomposePar
mpirun --hostfile machines -np 4 advectionFoam -parallel -heun2
reconstructPar
setAnalyticTracerField -time 5000,10000
#gnuplot -p stencil.plt

