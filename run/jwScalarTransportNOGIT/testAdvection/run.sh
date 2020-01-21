#!/bin/bash -e

rm -rf constant/polyMesh 0

blockMesh
#setFields
mkdir 0
setAnalyticTracerField
setVelocityField
cp 0/T_analytic 0/T

#scalarTransportFoamCN 
jwScalarTransportFoam -CN
#paraFoam

#-help can be used to see the options
