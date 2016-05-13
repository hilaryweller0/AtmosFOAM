#!/bin/bash -ve

# Create the case, run and post-process

rm -rf constant/polyMesh [0-9]* log gmt.history *.dat
blockMesh
setPlaneJet

# change zeroGradient boundary conditions of h to geostrophicZeroFlux to 
# ensure no flow over north and south boundaries despite a pressure graient
# in geostrophic balance with the wind
sed -i 's/zeroGradient/geostrophicZeroFlux/g' 0/h

gmtFoam -time 0 hU
evince 0/hU.pdf &

# Solve the shallow water equations
shallowWaterFoamBeta_C_explicit |& tee log

# Plot the final solutions
time=1e+06
gmtFoam -time $time hU
evince $time/hU.pdf &

# Plot animation of solutions
gmtFoam hU
animate 0/hU.pdf [1-9]00000/hU.pdf 1e+06/hU.pdf &

# Plot change between initial and final solutions
time=1e+06
sumFields $time hDiff $time h 0 h -scale0 1 -scale1 -1
sumFields $time UDiff $time U 0 U -scale0 1 -scale1 -1
gmtFoam -time $time hUDiff
evince $time/hUDiff.pdf &

# If we have an analytic solution in directory ref, then we can calculate 
# error norms relative to this solution. For demonstration, I assume that
# the analytic solution is the initial condition. This is not the case here
# as this is not a steady state problem
ref=0
time=1e+06
globalSum h -time $ref
globalSum hDiff -time $time
l1=`tail -1 globalSumh.dat | awk '{print $2}'`
l2=`tail -1 globalSumh.dat | awk '{print $3}'`
li=`tail -1 globalSumh.dat | awk '{print $4}'`
echo "time l1 l2 linf" > hErrors.dat
awk '{if ($1 == 1*$1) print $1, $2/'$l1', $3/'$l2', $4/'$li'}' \
    globalSumhDiff.dat >> hErrors.dat

# Calculate and plot the vorticity
time=1e+06
vorticity -time $time
writeuvw vorticity -time $time
gmtFoam -time $time vorticity
evince $time/vorticity.pdf &

