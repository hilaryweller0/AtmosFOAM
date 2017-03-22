#!/bin/bash
set -e
rm -rf [0-9]* processor*
setInitialTracerField
setVelocityField
decomposePar
mpirun --hostfile machines -np 2 advectionFoam -parallel -heun2
reconstructPar
setAnalyticTracerField -time 5000,10000
