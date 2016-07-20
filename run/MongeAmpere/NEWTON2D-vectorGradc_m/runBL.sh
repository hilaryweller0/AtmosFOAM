#!/bin/bash -e

rm -rf core constant [0-9]* log gmt.history convBL.eps
blockMesh
mkdir -p constant/rMesh 0
ln -s ../../gmtDicts constant/gmtDicts
cp init0/Phi 0
cp -r constant/polyMesh constant/rMesh
touch log
tail -f log &
BL2D >& log
pkill tail 

#exit

#gmtFoam -latestTime -region rMesh mesh
#gv [1-9]*/mesh.pdf &
#gmtFoam -latestTime Phi
#gv [1-9]*/Phi.pdf &

#gmtFoam -latestTime -region rMesh mesh
#gv [1-9]*/mesh.pdf &
grep PABe log | awk '{print $3, $6}' \
     | psxy -JX10c/7cl -R0/200/1e-8/10 -A \
         -Ba20:"Iterations":/a20:"Equidistribution": -W | ps2eps > convBL.eps
gv convBL.eps &

exit

time=5
gmtFoam -time $time lapc_m
gv $time/lapc_m.pdf &
gmtFoam -time $time Phi
gv $time/Phi.pdf &

