#!/bin/bash
set -e

rm -rf [0-9]*
blockMesh
mkdir -p 0
cp init_0/Uf 0/
setTheta
setExnerBalancedH
exnerFoamH

for t in [0-9]*
do
	sumFields -scale0 1 -scale1 -1 $t theta_diff $t theta 0 theta
done
