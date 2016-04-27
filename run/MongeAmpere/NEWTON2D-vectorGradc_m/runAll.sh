#!/bin/bash -e

# setup computational mesn and initial conditions mesh
rm -rf core constant [0-9]* log gmt.history conv.eps
mkdir -p save

# create a special case
case=save/60x60_bell_NewtonLinear_gradcmLinear
case=save/60x60_bell_NewtonLinear_gradcmDown
case=save/60x60_bell_NewtonUp_gradcmDown
case=save/100x100_bell_NewtonLinear_gradcmDown
case=save/200x200_bell_NewtonLinear_gradcmDown
case=save/200x200_bell_NewtonLinear_gradcmDown_relax6
case=save/200x200_bell_NewtonLinear_gradcmDown_relax2
case=save/300x300_bell_NewtonLinear_gradcmDown
case=save/300x300_bell_AFP
case=save/60x60_ring_NewtonLinear_gradcmDown
case=save/100x100_ring_NewtonLinear_gradcmDown
case=save/200x200_ring_NewtonLinear_gradcmDown
case=save/300x300_ring_NewtonLinear_gradcmDown
case=save/60x60_bell_NewtonLinear_gradcmDown_setRef
mkdir $case
cp -r system $case
gedit -s $case/system/fv* $case/system/*Dict  $case/system/rMesh/fvSchemes

# create the mesh and initial conditions
blockMesh -case $case
mkdir -p $case/constant/rMesh $case/0
ln -s ../../../../gmtDicts $case/constant/gmtDicts
cp init0/Phi $case/0
cp -r $case/constant/polyMesh $case/constant/rMesh

# Solve the MA equation
NEWTON2D-vectorGradc_m -case $case >& $case/log &
tail -f $case/log

# Diagnostics

grep PABe $case/log | awk '{print $3, $6}' \
     | psxy -JX10c/7cl -R0/200/1e-8/10 -A \
         -Ba20:"Iterations":/a20:"Equidistribution": -W | \
          ps2eps > $case/conv.eps
gv $case/conv.eps &
grep PABe $case/log | awk '{print $3, $6}' \
     | psxy -JX8.4c/6.3cl -R0/450/1e-8/10 -A -W \
         -Bf450/f10 -W \
         | ps2eps > $case/conv2.eps
gv $case/conv2.eps &

gmtFoam -case $case -latestTime -region rMesh mesh
gv $case/[1-9]*/mesh.pdf &
gmtFoam -case $case -latestTime Phi
gv $case/[1-9]*/Phi.pdf &
gmtFoam -case $case -latestTime lapc_m
gv $case/[1-9]*/lapc_m.pdf &

