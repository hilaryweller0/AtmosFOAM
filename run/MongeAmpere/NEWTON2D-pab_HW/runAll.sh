#!/bin/bash -e

rm -rf core constant [0-9]* log conv.eps gmt.history
blockMesh
mkdir -p constant/rMesh 0
ln -s ../../gmtDicts constant/gmtDicts
cp init0/Phi 0
cp -r constant/polyMesh constant/rMesh
touch log
tail -f log &
NEWTON2D-pab >& log
pkill tail 

#gmtFoam -latestTime -region rMesh mesh
#gv [1-9]*/mesh.pdf &
grep PABe log | awk '{print $3, $6}' \
     | psxy -JX10c/7cl -R0/200/1e-8/10 -A \
         -Ba20:"Iterations":/a20:"Equidistribution": -W | ps2eps > conv.eps
gv conv.eps &

