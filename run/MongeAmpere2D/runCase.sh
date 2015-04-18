#!/bin/bash -e

# generate the mesh
blockMesh
# copy the mesh to the re-distribued mesh
mkdir -p constant/rMesh
cp -r constant/polyMesh constant/rMesh

# solve Monge Ampere, pipting the output to a log file
MongeAmpere2D | tee log

# The numerical solution, Phi, on the old mesh is now in the time
# directories, the points of the new mesh for each itereation is in
# rMesh/polyMesh for each time step

# extract convergence diagnostics from the log file
echo "#time minSource maxSource minVol maxVol ratio boostLaplacian" > MAconvergence.dat
grep boostLaplacian log | awk '{print $3, $7, $9, $14, $16, $16/$14, $19}'\
    >> MAconvergence.dat
echo "resid nIters" > residsTmp.dat
grep residual log | awk '{print $8, $15}' | awk -F',' '{print $1, $2}'\
     >> residsTmp.dat
paste MAconvergence.dat residsTmp.dat > convAll.dat
rm residsTmp.dat
mv convAll.dat MAconvergence.dat

# plot some convergence diagnostics
gmtPlot plots/MAconvergence.gmt
gmtPlot plots/MAratios.gmt
gmtPlot plots/boostLaplacian.gmt
gmtPlot plots/MAresids.gmt
gmtPlot plots/nIterations.gmt

# find the final time (for the new mesh)
export time=`ls $case | sort -n | tail -1`

# Plot the solution of the MA equation and its gradient on the old mesh
gmtFoam -time $time Phi
evince $time/Phi.pdf &

# Plot the monitor function on the new mesh
gmtFoam -region rMesh -time $time monitor
evince $time/monitor.pdf &

# mesh analysis
export time=`ls $case | sort -n | tail -1`
meshAnalysis2D MongeAmpere2DDict -time $time -region rMesh
gmtPlot plots/distDx.gmt
gmtPlot plots/distArea.gmt

